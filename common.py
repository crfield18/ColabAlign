'''Shared constants and PDB parsing utilities, used by both structure_prep.py
and consensus.py so they can't silently drift apart from each other.'''

AA_STANDARD = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}

AA_MODIFIED_USALIGN_COMPATIBLE = {'MSE'}

AA_CODES = AA_STANDARD | AA_MODIFIED_USALIGN_COMPATIBLE

# selenomethionine (MSE) is compatible with US-align using the flag: -het 2
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M',
}

def chunks(lst, n):
    '''Yield n roughly-equal chunks from the list, for distributing work
    across a fixed number of worker processes.'''
    avg = len(lst) / n
    last = 0
    while last < len(lst):
        yield lst[int(last):int(last + avg)]
        last += avg


def file_contains_protein(pdb_path):
    '''Check both ATOM and HETATM lines, since USalign-alignable modified
    residues like MSE are conventionally recorded as HETATM.'''
    with open(pdb_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[17:20].strip() in AA_CODES:
                return True
    return False


def extract_sequence(pdb_path):
    '''Extract one-letter sequence from CA atoms of the first chain encountered.'''
    seq, seen_chain, seen_resnum = [], None, set()
    with open(pdb_path, 'r', encoding='utf8') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if line[12:16].strip() != 'CA':
                continue
            resname = line[17:20].strip()
            chain_id = line[21].strip()
            resnum = line[22:27].strip()
            if seen_chain is None:
                seen_chain = chain_id
            if chain_id != seen_chain:
                break
            if resnum in seen_resnum:
                continue
            seen_resnum.add(resnum)
            seq.append(THREE_TO_ONE.get(resname, 'X'))
    return ''.join(seq)
