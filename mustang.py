from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
from collections import defaultdict
from bs4 import BeautifulSoup
from Bio import SeqIO
from math import inf

class GetRepresentatives():
    def __init__(self, clusters: dict, threads: int = 1) -> None:
        self.threads = threads
        self.clusters = clusters

    def _run_mustang(self, pdb_dir: Path, pdb_names: list, cluster_num: int | str, threshold: int | str):
        cmd = [
            'mustang',
            '-p', str(pdb_dir),
            '-o', pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}_mustang'),
            '-F',  'fasta',
            '-s', 'OFF'
        ]

        cmd.extend(['-i'])
        cmd.extend(pdb_names)  # This will add each PDB name as a separate argument so MUSTANG won't throw a fit

        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            try:
                with open(pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}_mustang.txt'), 'w') as mustang_log:
                    stdout, stderr = process.communicate()
                    mustang_log.write(stdout.decode())
            finally:
                process.terminate()
        return stdout, stderr

    def _run_mview(self, input_alignment: Path):
        # https://desmid.github.io/mview/manual/manual.html
        cmd = [
            'mview',
            '-in', 'fasta',
            '-html', 'head',
            '-consensus', 'on',
            '-bold',
            '-css', 'on',
            '-coloring', 'any',
            '-con_coloring', 'any',
            '-con_threshold', '100,90,80,70,60',
            str(input_alignment)
        ]

        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            try:
                with open(input_alignment.with_suffix('.html'), 'w') as mview_file:
                    stdout, stderr = process.communicate()
                    mview_file.write(stdout.decode())
            finally:
                process.terminate()
        return stdout, stderr
    
    def _parse_mview(self, html_file: Path):
        assert html_file.suffix == '.html', 'Input file must be an HTML file'
        assert html_file.exists(), f'File does not exist: {html_file}'
        
        with open(html_file, 'r', encoding='utf-8') as in_file:
            html_content = in_file.read()

        soup = BeautifulSoup(html_content, 'html.parser')
        lines = soup.get_text().split('\n')

        with open(html_file.with_name(f'{html_file.stem}_consensus.afasta'), 'w', encoding='utf-8') as out_file:
            for line in lines:
                if not line.lstrip().startswith('consensus'):
                    continue

                name, sequence = line.strip().split(maxsplit=1)
                out_file.write(f'>{name}\n')
                for i in range(0, len(sequence), 80):
                    out_file.write(f'{sequence[i:i+80]}\n')

    def _find_nearest_sequence(self, mustang_afasta: Path, mview_afasta: Path):
        def _is_match(c1, c2):
            match (c1, c2):
                # Both characters are gap characters
                case ('.', '-'):
                    return True
                # If there is a gap in the consensus and no character in the query sequence
                case ('.', _):
                    return False
                # If the characters are identical
                case (x, y) if x == y:
                    return True
                # If the query character is in the mview look-up-table group for the consensus character
                case (x, y) if y in mview_LUT.get(x, set()):
                    return True
                # All other cases
                case _:
                    return False

        # Possible amino acid groupings in the MVIEW consensus sequence  output
        mview_LUT = {
        'o': {'S', 'T'},                                                              # alcohol
        'l': {'I', 'L', 'V'},                                                         # aliphatic
        'a': {'F', 'H', 'W', 'Y'},                                                    # aromatic
        'c': {'D', 'E', 'H', 'K', 'R'},                                               # charged
        'h': {'A', 'C', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'R', 'T', 'V', 'W', 'Y'},  # hydrophobic
        '-': {'D', 'E'},                                                              # negative
        'p': {'C', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'},                      # polar
        '+': {'H', 'K', 'R'},                                                         # positive
        's': {'A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'},                           # small
        'u': {'A', 'G', 'S'},                                                         # tiny
        't': {'A', 'C', 'D', 'E', 'G', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'},            # turnlike
        '*': {'*'}                                                                    # stop
        }

        closest_seq = {
            'record_name': '',
            'sequence': '',
            'score': -inf,
            'max_score': -inf
            }

        with open(mview_afasta, 'r', encoding='UTF8') as consensus_fasta:
            # Using the 60% consensus sequences to have the largest number of non-gap
            # residues as possible
            consensus = str(list(SeqIO.parse(consensus_fasta, 'fasta'))[-1].seq)
        
        with open(mustang_afasta, 'r', encoding='UTF8') as cluster_aligned_fasta:
            records = list(SeqIO.parse(cluster_aligned_fasta, 'fasta'))
            for record in records:
                score = 0
                # Need to change - to * to match the gap character in the BLOSUM matrix
                for i, (cons_char, seq_char) in enumerate(zip(consensus, record.seq)):
                    if _is_match(cons_char, seq_char):
                        score += 1
                    else:
                        score -= 1
                
                if score > closest_seq.get('score') or closest_seq.get('score') is -inf:
                    closest_seq = {
                        'record_name': record.name,
                        'sequence': record.seq,
                        'score': score,
                        'max_score': len(consensus)
                        }

        return closest_seq

    def _run_mustang_and_mview(self, pdb_dir: Path, pdb_names: list, cluster_num: int | str, threshold: int | str):
        print(f'Calculating representative for cluster:\t{cluster_num}')
        if not pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}_mustang.afasta').exists():
            stdout, stderr = self._run_mustang(pdb_dir, pdb_names, cluster_num, threshold)
            if stderr:
                print(f'Error in MUSTANG: {stderr.decode()}')
                return stdout, stderr

        alignment_file = pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}_mustang.afasta')
        stdout, stderr = self._run_mview(alignment_file)
        if stderr:
            print(f'Error in MView: {stderr.decode()}')
        
        self._parse_mview(pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}_mustang.html'))
        nearest_sequence_to_consensus = self._find_nearest_sequence(
            mustang_afasta=alignment_file,
            mview_afasta=alignment_file.with_name(f'{alignment_file.stem}_consensus.afasta')
            )

        return stdout, stderr, nearest_sequence_to_consensus, cluster_num, threshold

    def get_representatives(self, results_path: Path):
        all_clusters = defaultdict(lambda: {
            'threshold': None,
            'cluster': None,
            'dir': None
        })

        i = 0
        clusters_dir = results_path.joinpath('clusters')

        # Collect all cluster directories
        for threshold_dir in [t for t in clusters_dir.iterdir() if t.is_dir()]:
            for cluster_dir in [dir for dir in threshold_dir.iterdir() if dir.is_dir() and dir.stem != '-1']:
                all_clusters[i]['threshold'] = threshold_dir.name
                all_clusters[i]['cluster'] = cluster_dir.stem
                all_clusters[i]['dir'] = cluster_dir
                i += 1

        representatives = defaultdict(
            lambda : {
            -inf: {
                'record_name': None,
                'sequence': None,
                'score': -inf,
                'max_score': -inf
            }
            }
            )

        # Run MUSTANG in parallel for each cluster
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(
                    self._run_mustang_and_mview,
                    pdb_dir=info['dir'],
                    pdb_names=[f'/{f.name}' for f in info['dir'].iterdir() if f.suffix == '.pdb'],
                    cluster_num=info['cluster'],
                    threshold=info['threshold']
                )
                for _, info in all_clusters.items()
            ]

            for future in as_completed(futures):
                _, stderr, cluster_representative, cluster, threshold = future.result()
                if stderr:
                    print(f'Error: {stderr.decode()}')
                representatives[int(cluster)] = cluster_representative

        with open(results_path.joinpath(f'clusters/{threshold}/threshold_{threshold}_cluster_representatives.fasta'), 'w', encoding='UTF8') as representative_fasta:
            for num, entry in dict(sorted(representatives.items(), key=lambda x: int(x[0]))).items():
                representative_fasta.write(f'> Cluster {num}: {entry.get('record_name')} (alignment score {entry.get('score')}/{entry.get('max_score')})\n{str(entry.get('sequence')).replace("-", "")}\n')
