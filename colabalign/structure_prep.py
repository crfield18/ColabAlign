'''Prepare input structures (chain splitting, .cif conversion, filtering) for alignment.'''

import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.PDB import PDBParser, PDBIO
from tqdm.auto import tqdm as tqdm_auto

from .common import file_contains_protein, chunks


STRUCTURE_SUFFIXES = ('.pdb', '.cif')


def _is_structure_file(path: Path) -> bool:
    '''True for .pdb/.cif files (case-insensitive), ignoring macOS ._ resource forks.'''
    return (path.is_file() and path.suffix.lower() in STRUCTURE_SUFFIXES
            and not path.name.startswith('._'))


def _is_junk_zip_member(name: str) -> bool:
    '''Skip macOS metadata (__MACOSX/, ._*) and hidden files bundled into zips.'''
    parts = Path(name).parts
    return any(part == '__MACOSX' or part.startswith('.') for part in parts)


def extract_zip_archive(zip_path: Path, extract_dir: Path):
    '''
    Extract only the .pdb/.cif members of a zip archive into
    extract_dir/<zip stem>/, preserving the archive's internal folder layout
    so identically named files in different folders don't overwrite each
    other. Returns the list of extracted structure files.
    '''
    target = extract_dir.joinpath(zip_path.stem)
    counter = 1
    while target.exists():  # two zips with the same name from different folders
        target = extract_dir.joinpath(f'{zip_path.stem}_{counter}')
        counter += 1
    target.mkdir(parents=True)

    extracted = []
    try:
        with zipfile.ZipFile(zip_path) as archive:
            for member in archive.infolist():
                if member.is_dir() or _is_junk_zip_member(member.filename):
                    continue
                if Path(member.filename).suffix.lower() not in STRUCTURE_SUFFIXES:
                    continue
                # ZipFile.extract sanitises absolute paths and '..' components
                extracted.append(Path(archive.extract(member, target)))
    except zipfile.BadZipFile as e:
        print(f'Skipping unreadable zip file {zip_path}: {e}')
        return []

    print(f'Extracted {len(extracted)} structure file(s) from {zip_path.name}.')
    return extracted


def discover_input_models(input_paths, extract_dir: Path):
    '''
    Collect .pdb/.cif files from a mix of files, directories and .zip
    archives (given directly or found inside an input directory). Zip
    contents are extracted to extract_dir, which the caller can delete once
    structure preparation is complete.
    '''
    model_list, zip_files = [], []
    for user_input in input_paths:
        if user_input.is_file():
            if user_input.suffix.lower() == '.zip':
                zip_files.append(user_input)
            elif _is_structure_file(user_input):
                model_list.append(user_input)
        elif user_input.is_dir():
            for f in user_input.iterdir():
                if f.is_file() and f.suffix.lower() == '.zip':
                    zip_files.append(f)
                elif _is_structure_file(f):
                    model_list.append(f)

    # Output chain files are named <stem><chain>.pdb, so two inputs sharing a
    # stem would collide (and 'first' mode silently drops the _dup copy).
    # Zips often contain identically named files in different folders
    # (e.g. ranked_0.pdb per AlphaFold run), so extracted files with a
    # clashing stem are renamed with their parent folder as a prefix.
    seen_stems = {m.stem.lower() for m in model_list}
    for zip_path in zip_files:
        for f in extract_zip_archive(zip_path, extract_dir):
            if not _is_structure_file(f):
                continue
            if f.stem.lower() in seen_stems:
                f = _rename_to_unique_stem(f, seen_stems)
            seen_stems.add(f.stem.lower())
            model_list.append(f)

    return sorted(model_list)


def _rename_to_unique_stem(path: Path, seen_stems) -> Path:
    '''Rename an extracted file to <parent>_<stem>[_N]<suffix> so its stem is unique.'''
    base = f'{path.parent.name}_{path.stem}'
    candidate, counter = base, 1
    while candidate.lower() in seen_stems:
        candidate = f'{base}_{counter}'
        counter += 1
    new_path = path.with_name(f'{candidate}{path.suffix}')
    print(f'Renamed {path.name} from zip to {new_path.name} to avoid a filename clash.')
    return path.rename(new_path)


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


def _safe_copy_destination(dst: Path, src_size: int) -> Path:
    '''
    Return a non-clobbering destination path if a file with the same name
    already exists but has a different size (i.e. is a genuinely different
    file, not a re-run of the same input). Same-size existing files are
    assumed identical and left alone.
    '''
    if not dst.exists():
        return dst
    if dst.stat().st_size == src_size:
        return dst  # treat as already-copied, no action needed

    stem, suffix, parent = dst.stem, dst.suffix, dst.parent
    counter = 1
    while True:
        candidate = parent.joinpath(f'{stem}_dup{counter}{suffix}')
        if not candidate.exists() or candidate.stat().st_size == src_size:
            return candidate
        counter += 1


def _safe_move(src: Path, dst: Path):
    '''Move src to dst, redirecting to a non-clobbering path if dst already exists.'''
    dst = _safe_copy_destination(dst, src.stat().st_size)
    if dst != src:
        shutil.move(str(src), str(dst))


def split_chains_to_pdb(models_path: Path, model: Path):
    '''Split a native PDB file into one file per chain using biopython.'''
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

        final_path = _safe_copy_destination(models_path.joinpath(out_name), tmp_path.stat().st_size)
        tmp_path.rename(final_path)


def _prepare_one_model(beem_path: Path, models_path: Path, model: Path):
    '''
    Convert/copy a single input model into single-chain PDB file(s) under
    models_path. Runs inside a worker process.
    Never raises - returns (model_stem, error_or_None) so that one bad
    input can't take down the rest of its chunk.
    '''
    try:
        if model.suffix.lower() == '.cif':
            # Private working dir per call: BeEM writes into its cwd, and
            # parallel workers must not see (or move) each other's outputs.
            with tempfile.TemporaryDirectory(dir=models_path) as work_dir:
                work_dir = Path(work_dir)
                stdout, stderr = _run_beem(beem_path, model, work_dir)
                beem_models = [item for item in stdout.decode(errors='ignore').split('\n') if item]
                if not beem_models:
                    return model.stem, f'BeEM produced no output (stderr: {stderr.decode(errors="ignore")})'
                for m in beem_models:
                    _safe_move(work_dir.joinpath(m), models_path.joinpath(Path(m).name))
        else:
            split_chains_to_pdb(models_path, model)
    except Exception as e:  # pylint: disable=broad-exception-caught
        # any failure (subprocess, parsing, etc.) is reported as a failed
        # model rather than crashing the whole chunk
        return model.stem, f'{type(e).__name__}: {e}'
    return model.stem, None


def _prepare_job(beem_path: Path, models_path: Path, models_chunk):
    '''Prepare every model in a chunk, sequentially, within one worker process.'''
    return [_prepare_one_model(beem_path, models_path, model) for model in models_chunk]


def _copy_and_convert_models(model_list, models_path: Path, beem_path: Path, cores: int):
    '''
    Get the full list of input models up front, split it into per-core
    chunks, and run BeEM conversion / chain splitting for each chunk in
    a separate worker process.
    '''
    chunk_count = min(int(cores), len(model_list)) or 1
    jobs = list(chunks(model_list, chunk_count))

    total_bar = tqdm_auto(total=len(model_list), desc='Preparing structures',
                           unit=' structure', dynamic_ncols=True, leave=True)

    failed = []
    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [executor.submit(_prepare_job, beem_path, models_path, job) for job in jobs]

        for future in as_completed(futures):
            try:
                results = future.result()
            except Exception as e:  # pylint: disable=broad-exception-caught
                # one lost chunk shouldn't abort collection of the rest
                print(f'Chunk-level error (entire chunk lost): {e}')
                continue

            for model_stem, error in results:
                total_bar.update(1)
                if error is not None:
                    failed.append((model_stem, error))

    total_bar.close()

    if failed:
        print(f'\n{len(failed)} input structure(s) failed during preparation and were skipped:')
        for stem, err in failed:
            print(f'  {stem}: {err}')


def prepare_models(model_list, models_path: Path, beem_path: Path, mode: str, cores: int = 1):
    '''
    Copy/convert input structures into single-chain PDB files under models_path.
    Returns the final, filtered list of PDB files to align.
    '''
    print('Copying pdb files and converting .cif files to .pdb.')
    original_stems = {m.stem for m in model_list}

    _copy_and_convert_models(model_list, models_path, beem_path, cores)

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

    return final_model_list
