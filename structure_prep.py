'''Prepare input structures (chain splitting, .cif conversion, filtering) for alignment.'''

import shutil
import subprocess
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO

from common import file_contains_protein


def discover_input_models(input_paths):
    '''Collect .pdb/.cif files from a mix of files and directories.'''
    model_list = []
    for user_input in input_paths:
        if user_input.is_file() and user_input.suffix in ('.pdb', '.cif'):
            model_list.append(user_input)
        elif user_input.is_dir():
            model_list.extend(
                f for f in user_input.iterdir()
                if f.is_file() and f.suffix in ('.pdb', '.cif')
            )
    return sorted(model_list)


def _run_beem(beem_path: Path, models_path: Path, cif_file: Path):
    '''Run BeEM to convert a .cif file into per-chain .pdb files.'''
    cmd = [
        str(beem_path),
        models_path.joinpath(cif_file.name).as_posix(),
        '-outfmt=2',
        f'-p={cif_file.stem}'
    ]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
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


def prepare_models(model_list, models_path: Path, beem_path: Path, mode: str):
    '''
    Copy/convert input structures into single-chain PDB files under models_path.
    Returns the final, filtered list of PDB files to align.
    '''
    print('Copying pdb files and converting .cif files to .pdb.')
    original_stems = {m.stem for m in model_list}

    for model in model_list:
        if model.suffix == '.cif':
            stdout, _ = _run_beem(beem_path, models_path, model)
            beem_models = (item for item in stdout.decode(errors='ignore').split('\n') if item)
            for m in beem_models:
                _safe_move(Path(m), models_path.joinpath(Path(m).name))
        else:
            split_chains_to_pdb(models_path, model)

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
