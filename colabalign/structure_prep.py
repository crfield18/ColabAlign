'''Prepare input structures (chain splitting, .cif conversion, filtering) for alignment.'''

import csv
import hashlib
import shutil
import subprocess
import tarfile
import tempfile
import zipfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.PDB import PDBParser, PDBIO
from tqdm.auto import tqdm as tqdm_auto

from .common import file_contains_protein, chunks, process_context


STRUCTURE_SUFFIXES = ('.pdb', '.cif')


def _is_structure_file(path: Path) -> bool:
    '''True for .pdb/.cif files (case-insensitive), ignoring macOS ._ resource forks.'''
    return (path.is_file() and path.suffix.lower() in STRUCTURE_SUFFIXES
            and not path.name.startswith('._'))


# Longest first so 'x.tar.gz' is matched as '.tar.gz', not '.gz'.
TAR_SUFFIXES = ('.tar.gz', '.tar.bz2', '.tar.xz', '.tgz', '.tbz2', '.txz', '.tar')
ARCHIVE_SUFFIXES = ('.zip',) + TAR_SUFFIXES


def _archive_suffix(path: Path):
    '''Return the (lower-case) archive suffix of path, or None if it is not an archive.'''
    name = path.name.lower()
    return next((sfx for sfx in ARCHIVE_SUFFIXES if name.endswith(sfx)), None)


def _is_archive(path: Path) -> bool:
    return path.is_file() and _archive_suffix(path) is not None


def _archive_stem(path: Path) -> str:
    ''''models.tar.gz' -> 'models' (Path.stem would give 'models.tar').'''
    return path.name[:len(path.name) - len(_archive_suffix(path))]


def _is_junk_member(name: str) -> bool:
    '''Skip macOS metadata (__MACOSX/, ._*) and hidden files bundled into archives.'''
    parts = Path(name).parts
    return any(part == '__MACOSX' or (part.startswith('.') and part not in ('.', '..'))
               for part in parts)


def _wanted_member(name: str) -> bool:
    return not _is_junk_member(name) and Path(name).suffix.lower() in STRUCTURE_SUFFIXES


def _unique_target(extract_dir: Path, stem: str) -> Path:
    target = extract_dir.joinpath(stem)
    counter = 1
    while target.exists():  # two archives with the same name from different folders
        target = extract_dir.joinpath(f'{stem}_{counter}')
        counter += 1
    target.mkdir(parents=True)
    return target


def _extract_zip_members(archive_path: Path, target: Path):
    '''Extract wanted zip members; returns [(extracted_path, member_name), ...].'''
    extracted = []
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            if member.is_dir() or not _wanted_member(member.filename):
                continue
            # ZipFile.extract sanitises absolute paths and '..' components
            extracted.append((Path(archive.extract(member, target)), member.filename))
    return extracted


def _extract_tar_members(archive_path: Path, target: Path):
    '''Extract wanted tar members; returns [(extracted_path, member_name), ...].'''
    extracted = []
    # 'r:*' auto-detects gzip/bzip2/xz/uncompressed.
    with tarfile.open(archive_path, 'r:*') as archive:
        for member in archive:
            # Regular files only: no symlinks, hardlinks or device files.
            if not member.isfile() or not _wanted_member(member.name):
                continue
            # The 'data' filter rejects absolute paths and '..' escapes
            # (raises tarfile.FilterError) and strips unsafe permissions.
            try:
                archive.extract(member, target, filter='data')
            except tarfile.FilterError as e:
                print(f'Skipping unsafe archive member {member.name}: {e}')
                continue
            extracted.append((target.joinpath(member.name), member.name))
    return extracted


def _extract_archive_members(archive_path: Path, extract_dir: Path):
    target = _unique_target(extract_dir, _archive_stem(archive_path))
    try:
        if _archive_suffix(archive_path) == '.zip':
            extracted = _extract_zip_members(archive_path, target)
        else:
            extracted = _extract_tar_members(archive_path, target)
    except (zipfile.BadZipFile, tarfile.TarError, EOFError, OSError) as e:
        print(f'Skipping unreadable archive {archive_path}: {e}')
        return []
    print(f'Extracted {len(extracted)} structure file(s) from {archive_path.name}.')
    return extracted


def content_hash(path: Path) -> str:
    '''Full sha256 hex digest of a file's contents.'''
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


DUPLICATES_CSV_COLUMNS = (
    'input_file_name', 'input_file_path', 'name_in_colabalign', 'sha256', 'identical_to',
)


def _iter_candidates(input_paths, extract_dir: Path):
    '''
    Yield (working_path, original_name, original_path, is_colabalign_copy)
    for every structure file: loose files first (in input order, directories
    sorted), then the contents of archives in the order they were found.
    Archive members are reported as '<archive path>:<member path>'.
    '''
    loose, archives, seen_paths = [], [], set()

    def add(path: Path, bucket):
        resolved = path.resolve()
        if resolved in seen_paths:  # the same file given twice, e.g. -i dir dir/a.pdb
            return
        seen_paths.add(resolved)
        bucket.append(path)

    for user_input in input_paths:
        user_input = Path(user_input)
        if user_input.is_file():
            if _is_archive(user_input):
                add(user_input, archives)
            elif _is_structure_file(user_input):
                add(user_input, loose)
        elif user_input.is_dir():
            for f in sorted(user_input.iterdir()):
                if _is_archive(f) and not f.name.startswith('._'):
                    add(f, archives)
                elif _is_structure_file(f):
                    add(f, loose)

    for path in loose:
        yield path, path.name, str(path.resolve()), False
    for archive_path in archives:
        for path, member in _extract_archive_members(archive_path, extract_dir):
            if _is_structure_file(path):
                yield path, Path(member).name, f'{archive_path.resolve()}:{member}', True


def collect_input_models(input_paths, extract_dir: Path):
    '''
    Collect .pdb/.cif files from a mix of files, directories and .zip/tar
    archives (given directly or found inside an input directory). Archive
    contents are extracted to extract_dir, which the caller can delete once
    structure preparation is complete.

    Output chain files are named <stem><chain>.pdb, so every input needs a
    unique stem (case-insensitive, so a.pdb and A.cif also clash). The first
    file with a given stem keeps its name; later ones get a working copy
    named <stem>_dup1, <stem>_dup2, ... Original files are never renamed:
    clashing loose files are copied into extract_dir, and clashing archive
    members are renamed in their extracted (temporary) copy.

    Returns (model_list, duplicate_rows). duplicate_rows has one dict per
    file in every clashing group (including the one that kept its name),
    keyed by DUPLICATES_CSV_COLUMNS.
    '''
    taken = set()     # lower-case stems in use
    groups = {}       # lower-case original stem -> [row, ...]
    model_list = []
    renamed_dir = None

    for path, orig_name, orig_path, is_copy in _iter_candidates(input_paths, extract_dir):
        stem = Path(orig_name).stem
        new_stem = _claim_unique_stem(stem, taken)
        if new_stem != stem:
            new_name = f'{new_stem}{path.suffix}'
            if is_copy:
                path = path.rename(path.with_name(new_name))
            else:
                if renamed_dir is None:
                    renamed_dir = _unique_target(extract_dir, '_renamed_inputs')
                path = Path(shutil.copy2(path, renamed_dir.joinpath(new_name)))
            print(f'Duplicate name: {orig_path} will be aligned as {new_name}.')

        model_list.append(path)
        groups.setdefault(stem.lower(), []).append({
            'input_file_name': orig_name,
            'input_file_path': orig_path,
            'name_in_colabalign': path.name,
            'sha256': content_hash(path),
            'identical_to': '',
        })

    duplicate_rows = []
    for rows in groups.values():
        if len(rows) > 1:
            _mark_identical(rows)
            duplicate_rows.extend(rows)
    return sorted(model_list), duplicate_rows


def _claim_unique_stem(stem: str, taken) -> str:
    '''Return stem, or stem_dup1, stem_dup2, ... (first one not in taken), and claim it.'''
    new_stem, counter = stem, 1
    while new_stem.lower() in taken:
        new_stem = f'{stem}_dup{counter}'
        counter += 1
    taken.add(new_stem.lower())
    return new_stem


def _mark_identical(rows):
    '''Fill identical_to with the first earlier row in the group that has the same sha256.'''
    first_with_hash = {}
    for row in rows:
        earlier = first_with_hash.setdefault(row['sha256'], row['name_in_colabalign'])
        if earlier != row['name_in_colabalign']:
            row['identical_to'] = earlier


def write_duplicates_csv(rows, csv_path: Path):
    '''Write the duplicate-name report (header only if there were no clashes).'''
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=DUPLICATES_CSV_COLUMNS, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def _run_beem(beem_path: Path, cif_file: Path, work_dir: Path):
    '''Run BeEM to convert a .cif file into per-chain .pdb files written to work_dir.'''
    cmd = [
        str(beem_path),
        cif_file.resolve().as_posix(),
        '-outfmt=2',
        f'-p={cif_file.stem}'
    ]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          cwd=work_dir) as process:
        stdout, stderr = process.communicate()
    return stdout, stderr


def _safe_copy_destination(dst: Path, src: Path) -> Path:
    '''
    Return a non-clobbering destination path if a file with the same name
    already exists with different contents (sha256). An existing file with
    identical contents is assumed to be a re-run of the same input and left
    alone.
    '''
    if not dst.exists():
        return dst
    src_hash = content_hash(src)
    if content_hash(dst) == src_hash:
        return dst  # treat as already-copied, no action needed

    stem, suffix, parent = dst.stem, dst.suffix, dst.parent
    counter = 1
    while True:
        candidate = parent.joinpath(f'{stem}_dup{counter}{suffix}')
        if not candidate.exists() or content_hash(candidate) == src_hash:
            return candidate
        counter += 1


def _safe_move(src: Path, dst: Path) -> Path:
    '''Move src to dst, redirecting to a non-clobbering path if dst already exists.
    Returns the path the file ended up at.'''
    dst = _safe_copy_destination(dst, src)
    if dst != src:
        shutil.move(str(src), str(dst))
    return dst


def split_chains_to_pdb(models_path: Path, model: Path):
    '''Split a native PDB file into one file per chain using biopython.
    Returns the names of the chain files written (or already in place).'''
    written = []
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure(model.stem, model)
    for chain in structure.get_chains():
        chain_id = chain.get_id().upper()
        out_name = f'{structure.get_id()}{chain_id}.pdb'
        io.set_structure(chain)

        # write to a temp path first so we can size-check before deciding
        # whether this collides with an existing file of the same name
        tmp_path = models_path.joinpath(f'.tmp_{out_name}')
        io.save(tmp_path.as_posix())

        final_path = _safe_copy_destination(models_path.joinpath(out_name), tmp_path)
        if final_path.exists():
            tmp_path.unlink()  # identical file already in place
        else:
            tmp_path.rename(final_path)
        written.append(final_path.name)
    return written


def _prepare_one_model(beem_path: Path, models_path: Path, model: Path):
    '''
    Convert/copy a single input model into single-chain PDB file(s) under
    models_path. Runs inside a worker process.
    Never raises - returns (model, error_or_None, [chain file names]) so
    that one bad input can't take down the rest of its chunk.
    '''
    outputs = []
    try:
        if model.suffix.lower() == '.cif':
            # Private working dir per call: BeEM writes into its cwd, and
            # parallel workers must not see (or move) each other's outputs.
            with tempfile.TemporaryDirectory(dir=models_path) as work_dir:
                work_dir = Path(work_dir)
                stdout, stderr = _run_beem(beem_path, model, work_dir)
                beem_models = [item for item in stdout.decode(errors='ignore').split('\n') if item]
                if not beem_models:
                    return model, f'BeEM produced no output (stderr: {stderr.decode(errors="ignore")})', []
                for m in beem_models:
                    dst = _safe_move(work_dir.joinpath(m), models_path.joinpath(Path(m).name))
                    outputs.append(dst.name)
        else:
            outputs = split_chains_to_pdb(models_path, model)
    except Exception as e:  # pylint: disable=broad-exception-caught
        # any failure (subprocess, parsing, etc.) is reported as a failed
        # model rather than crashing the whole chunk
        return model, f'{type(e).__name__}: {e}', outputs
    return model, None, outputs


def _prepare_job(beem_path: Path, models_path: Path, models_chunk):
    '''Prepare every model in a chunk, sequentially, within one worker process.'''
    return [_prepare_one_model(beem_path, models_path, model) for model in models_chunk]


def _copy_and_convert_models(model_list, models_path: Path, beem_path: Path, cores: int):
    '''
    Get the full list of input models up front, split it into per-core
    chunks, and run BeEM conversion / chain splitting for each chunk in
    a separate worker process.
    Returns {input model path: [chain file names it produced]}.
    '''
    chunk_count = min(int(cores), len(model_list)) or 1
    jobs = list(chunks(model_list, chunk_count))

    total_bar = tqdm_auto(total=len(model_list), desc='Preparing structures',
                           unit=' structure', dynamic_ncols=True, leave=True)

    failed = []
    outputs_by_model = {}
    with ProcessPoolExecutor(max_workers=cores, mp_context=process_context()) as executor:
        futures = [executor.submit(_prepare_job, beem_path, models_path, job) for job in jobs]

        for future in as_completed(futures):
            try:
                results = future.result()
            except Exception as e:  # pylint: disable=broad-exception-caught
                # one lost chunk shouldn't abort collection of the rest
                print(f'Chunk-level error (entire chunk lost): {e}')
                continue

            for model, error, outputs in results:
                total_bar.update(1)
                outputs_by_model[model] = outputs
                if error is not None:
                    failed.append((model.stem, error))

    total_bar.close()

    if failed:
        print(f'\n{len(failed)} input structure(s) failed during preparation and were skipped:')
        for stem, err in failed:
            print(f'  {stem}: {err}')
    return outputs_by_model


def prepare_models_with_chains(model_list, models_path: Path, beem_path: Path, mode: str,
                               cores: int = 1):
    '''
    Copy/convert input structures into single-chain PDB files under models_path.
    Returns (the final, filtered list of PDB files to align,
    {input file name: [aligned chain labels]}), where a label is the chain file
    stem used in the score matrix, tree and clusters (e.g. ranked_0_dup1A).
    Inputs with no surviving chains map to [].
    '''
    print('Copying pdb files and converting .cif files to .pdb.')
    original_stems = {m.stem for m in model_list}

    outputs_by_model = _copy_and_convert_models(model_list, models_path, beem_path, cores)

    # Drop empty or non-protein files - US-align does not handle them properly.
    for structure_file in models_path.glob('*.pdb'):
        if not structure_file.is_file() or structure_file.stat().st_size == 0:
            print(f'Removing empty file: {structure_file}')
            structure_file.unlink()
            continue

        if not file_contains_protein(structure_file):
            print(f'Removing non-protein file: {structure_file}')
            structure_file.unlink()

    final_model_list = sorted(models_path.glob('*.pdb'))

    if mode == 'first':
        seen_inputs = set()
        first_chains = []
        for model in final_model_list:
            input_stem = model.stem[:-1]
            if input_stem in original_stems and input_stem not in seen_inputs:
                seen_inputs.add(input_stem)
                first_chains.append(model)
            else:
                model.unlink()
        final_model_list = first_chains

    kept = {m.name for m in final_model_list}
    chains_by_input = {
        model.name: sorted(Path(n).stem for n in outputs_by_model.get(model, []) if n in kept)
        for model in model_list
    }
    return final_model_list, chains_by_input


def expand_duplicate_rows(rows, chains_by_input):
    '''
    Turn the per-file duplicate rows from collect_input_models() into one row
    per aligned chain, so name_in_colabalign matches the score matrix, tree
    and cluster labels (ranked_0_dup1A, ranked_0_dup1B, ...). identical_to
    points at the same chain of the earlier identical file. An input that
    produced no alignable chain keeps a single row with a blank
    name_in_colabalign.
    '''
    stem_of = {r['name_in_colabalign']: Path(r['name_in_colabalign']).stem for r in rows}
    expanded = []
    for row in rows:
        own_stem = stem_of[row['name_in_colabalign']]
        chains = chains_by_input.get(row['name_in_colabalign'], [])
        if not chains:
            expanded.append({**row, 'name_in_colabalign': '', 'identical_to': ''})
            continue
        earlier = row['identical_to']
        earlier_chains = set(chains_by_input.get(earlier, []))
        for label in chains:
            match = ''
            if earlier:
                candidate = stem_of[earlier] + label[len(own_stem):]
                match = candidate if candidate in earlier_chains else ''
            expanded.append({**row, 'name_in_colabalign': label, 'identical_to': match})
    return expanded
