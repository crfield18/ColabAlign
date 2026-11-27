# ColabAlign

## Fast pairwise protein secondary structure comparisons using multiprocessing

[![License](https://anaconda.org/bioconda/colabalign/badges/license.svg)](https://github.com/crfield18/ColabAlign/blob/main/LICENSE)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/colabalign/badges/platforms.svg)](https://anaconda.org/bioconda/colabalign)
[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/colabalign/README.html)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/crfield18/ColabAlign/blob/main/colabalign.ipynb)
[![Code DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.14169501-yellow)](https://doi.org/10.5281/zenodo.14169501)

This notebook performs pairwise protein structural alignments using the [US-align](https://zhanggroup.org/US-align/) algorithm by [Zhang _et al._, (2022)](https://doi.org/10.1038/s41592-022-01585-1), then constructs a structure-informed dendrogram using the UPGMA algorithm to visualise similarities.

ColabAlign is designed to run directly in Google Colab for ease-of-use and to remove any local hardware requirements. This implementation also includes multiprocessing support for dramatically increased performance over the base US-align program.

---
### Interpreting MView-like alignments

MView uses an expanded character set to represent groups of amino acids with similar properties. Below is 
Further information can be found at: <https://desmid.github.io/mview/manual/manual.html>

 **Group** | **MView character** | **1 letter amino acid codes**             
-------------------|------|---------------------------------------------------
 **Alcohol**       | `o`  | `S, T`                                     
 **Aliphatic**     | `l`  | `I, L, V`                                  
 **Aromatic**      | `a`  | `F, H, W, Y`                               
 **Charged**       | `c`  | `D, E, H, K, R`                            
 **Hydrophobic**   | `h`  | `A, C, F, G, H, I, K, L, M, R, T, V, W, Y` 
 **Negative**      | `-`  | `D, E`                                     
 **Polar**         | `p`  | `C, D, E, H, K, N, Q, R, S, T`             
 **Positive**      | `+`  | `H, K, R`                                  
 **Small**         | `s`  | `A, C, D, G, N, P, S, T, V`                
 **Tiny**          | `u`  | `A, G, S`                                  
 **Turn-like**     | `t`  | `A, C, D, E, G, H, K, N, Q, R, S, T`       
 **Stop**          | `*`  | `*`                                        


---

### Installation for local usage


**Conda package (once published).** Installs ColabAlign and all external tools:

`conda install -c conda-forge -c bioconda colabalign=1.2.3`

**From source.** Create the environment from the provided [YAML file](colabalign.yml), which includes the external tools, then install the package:

```
conda env create -f colabalign.yml
conda activate colabalign
python -m pip install --no-deps --no-build-isolation .
```

For exact reproduction including transitive dependencies and binary builds, use the platform-specific explicit locks instead of resolving `colabalign.yml` afresh. For example, on Linux x86-64:

```bash
conda create -n colabalign --file locks/conda-linux-64.lock
conda activate colabalign
python -m pip install --no-deps --no-build-isolation .
```

### Usage

```
colabalign -i structures/ more.zip runs.tar.gz extra.pdb -o results_dir -c 8 -n 50 -m first
```

`-n/--panel-size` is the number of structures you want to select. ColabAlign finds the clustering threshold that gives the closest achievable panel size. Not every size is possible, so the panel can be a little larger or smaller than requested. Give a range such as `-n 40-60` to get every achievable panel size in that range, or several values (`-n 30 50`). `results/panel_sizes.csv` lists every achievable size. You can still set thresholds directly with `-t 0.20 0.30`. With neither option, `-t 0.20` is used.

By default 5 % of each panel is given to **low-coverage structures**: the structures with the lowest best TM-score to the panel, added one at a time. They catch outliers inside large clusters that the cluster representative does not cover well. Change the share with `-l/--low-coverage-fraction` (for example `-l 0.1`), or turn it off with `-l 0`. Each panel is written to `results/clusters/<t>/threshold_<t>_panel.csv` and `.fasta`, and `results/panels.csv` summarises them all. Add `--no-superpose` to skip writing superposed structure files when you only need the panels; this is quicker, and the representatives and panels are the same. See the [manual](colabalign/MANUAL.md) for details.

`-i` accepts `.pdb`/`.cif` files, directories of them, and `.zip` or tar archives (`.tar`, `.tar.gz`/`.tgz`, `.tar.bz2`, `.tar.xz`). `python -m colabalign` works too. From Python:

```python
from colabalign import ColabAlign

if __name__ == '__main__':
    ColabAlign(inputs=['structures/'], output='results_dir', cores=8, panel_size=['40-60']).run()
```

#### Inputs with the same name

Every structure needs a unique name (ignoring case and extension, so `model.pdb` and `MODEL.cif` clash). When two inputs share a name, for example `ranked_0.pdb` from several AlphaFold runs, the first one found keeps its name and later ones are aligned as `ranked_0_dup1`, `ranked_0_dup2`, and so on. Loose files are taken first, in the order given to `-i` (directory contents in alphabetical order), then archive contents. Your original files are never renamed. ColabAlign works on renamed copies that are deleted once preparation finishes.

Each run writes `results/duplicates.csv`, which lists every file in a clashing group, including the one that kept its name. Files are split into chains before alignment, so there is one row per aligned chain, and `name_in_colabalign` matches the score matrix, tree and cluster labels:

| column | meaning |
|---|---|
| `input_file_name` | original file name |
| `input_file_path` | absolute path, or `<archive path>:<path inside archive>` for archived files |
| `name_in_colabalign` | chain label used in the run, e.g. `ranked_0_dup1A`, `ranked_0_dup1B`. Blank if the file produced no alignable protein chain |
| `sha256` | hash of the whole input file |
| `identical_to` | if this file is byte-identical to an earlier file in the group, the matching chain label from that file |

The CSV only has its header row when there were no clashes.

Keep the `if __name__ == '__main__':` guard in scripts. ColabAlign runs work in separate processes, and on every supported Python version (3.12 or newer, Linux and macOS) those processes import the calling script, so an unguarded call would start the run again inside each worker and fail. The command line, `python -m colabalign` and notebooks need no guard.

---

### Structure-informed MSA distance weights

The structure-informed MSA weights US-align residue correspondences by $w_k = \dfrac{1}{[1+(d_k/d_0)^2]}$, using the shorter full chain length reported from the supplied coordinates, rather than the aligned-pair count. Following the [US-align protein distance-scale calculation](https://raw.githubusercontent.com/pylelab/USalign/master/param_set.h), $d_0=0.5~\text{Å}$ through 21 residues and $d_0=1.24(L-15)^{1/3}-1.8~\text{Å}$ otherwise.

Individual weights remain between 0 and 1 and are not divided by chain length. This changes the MSA distance scale, not TM-score normalisation; the reported TM-score matrix, tree/clustering procedure and 6/1 profile gap penalties are unchanged. MSAs and selected representatives can change, without implying improved alignment accuracy.

Pair caches use a `v2_` filename namespace, so older weights are recomputed without deleting the old files. Missing or inconsistent distance tables retain the existing 1.0/0.5 match-symbol fallback; missing or nonpositive chain lengths fail the pair rather than silently using a different scale.

### Running the tests

```bash
python -m pip install --no-deps --no-build-isolation -e ".[test]"
pytest
```

The reference Conda environment already includes the pinned test tools. Tests that need USalign, BeEM or TreeCluster.py are skipped when those tools are not on `PATH`; release validation should run with all three installed.

### References

BibTeX-formatted references for this project and the associated references can be found in [references.bib](references.bib).

- **ColabAlign (paper)**: <https://doi.org/10.1101/2025.10.06.677802>
- **ColabAlign (code)**: <https://doi.org/10.5281/ZENODO.14169501>
- **BioPython**: <https://doi.org/10.1093/bioinformatics/btp163>
- **SciPy**: <https://doi.org/10.1038/s41592-019-0686-2>
- **US-align**: <https://doi.org/10.1038/s41592-022-01585-1>
- **TreeCluster**: <https://doi.org/10.1371/journal.pone.0221068>
- **MView**: <https://doi.org/10.1093/bioinformatics/14.4.380>
- **BeEM**: <https://doi.org/10.1186/s12859-023-05388-9>

---

> [!NOTE]
> Permission to use, copy, modify, and distribute this program for any purpose, with or without fee, is hereby granted, provided that the notices on the head, the reference information, and this copyright notice appear in all copies or substantial portions of the Software. It is provided "as is" without express or implied warranty.
