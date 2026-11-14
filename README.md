# ColabAlign

## Fast pairwise protein secondary structure comparisons using multiprocessing

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/crfield18/ColabAlign/blob/main/colabalign.ipynb) [![Code DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.14169501-blue)](https://doi.org/10.5281/zenodo.14169501) [![Paper DOI](https://img.shields.io/badge/DOI-10.1101/2025.10.06.677802-green)](https://doi.org/10.1101/2025.10.06.677802)

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

ColabAlign is a Python package with a `colabalign` command. It also needs three external tools: [USalign](https://github.com/pylelab/USalign), [BeEM](https://github.com/kad-ecoli/BeEM/) (only for `.cif` inputs) and [TreeCluster](https://github.com/niemasd/TreeCluster).

The v1.2.0 reference environment uses Python 3.12.16 and exact dependency versions. See [the dependency and release notes](DEPENDENCIES.md) for the pinning policy, complete environment locks and release checks. Biopython is pinned to 1.88, which includes the corrected UPGMA clade-height calculation; historical trees from older Biopython versions may therefore differ.

**Conda package (once published).** Installs ColabAlign and all external tools:

`conda install -c conda-forge -c bioconda colabalign=1.2.0`

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

**pip only.** `pip install .` installs the exactly pinned direct Python dependencies. Use the Python lock as described in [DEPENDENCIES.md](DEPENDENCIES.md) to pin transitive dependencies too. USalign 20241201 and BeEM 1.0.1 must be installed separately. Point ColabAlign at the binaries with `-u /path/to/USalign` and `-b /path/to/BeEM` if they are not on your `PATH`.

### Usage

```
colabalign -i structures/ more.zip extra.pdb -o results_dir -c 8 -t 0.20 0.30 -m first
```

`-i` accepts `.pdb`/`.cif` files, directories of them, and `.zip` archives. `python -m colabalign` works too. From Python:

```python
from colabalign import ColabAlign
ColabAlign(inputs=['structures/'], output='results_dir', cores=8, threshold=[0.2, 0.3]).run()
```

---

### Structure-informed MSA distance weights

The v1.2.0 structure-informed MSA uses US-align residue correspondences with distance weights \(w_k = 1/[1+(d_k/d_0)^2]\). The length used to calculate \(d_0\) is the shorter of the two **full chain lengths reported by US-align**, not the number of aligned residue pairs. For proteins, \(d_0=0.5\) Å for lengths up to 21 residues; otherwise \(d_0=1.24(L-15)^{1/3}-1.8\) Å, following the [US-align protein distance-scale calculation](https://raw.githubusercontent.com/pylelab/USalign/master/param_set.h).

This is a symmetric choice of distance scale for the MSA library, not a new implementation of a complete TM-score. The individual weights remain between 0 and 1 and are **not divided by chain length**. Profile gap-opening and gap-extension penalties remain 6 and 1. The two reported TM-scores, their retained maximum, the structural dendrogram and clustering procedure are unchanged.

Earlier implementations based \(d_0\) on the aligned-pair count. Changing to a full-chain scale can change MSAs, consensuses and selected representatives, particularly for partial alignments; it is an intentional algorithm revision, not a claim of improved alignment accuracy. Performance benchmarks and representative-selection comparisons should record the exact source commit and distinguish this revision from batching-only changes.

Per-pair caches now use a `v2_` filename namespace. Unversioned and older-version cached weights are ignored and recomputed automatically, without deleting the older files. Valid v2 results remain reusable across interrupted runs and different batching layouts. When distances cannot be parsed consistently, the existing match-symbol fallback weights of 1.0 and 0.5 are retained.

Relevant limits and edge cases are:

- **Partial alignments and long tails.** The spatial tolerance no longer shrinks when fewer residues align. The shorter-chain convention limits the effect of one much longer partner, but two long multidomain proteins or proteins with large unaligned tails can still receive a permissive distance scale. This is not an alignment-coverage filter.
- **Incomplete structures.** Full chain length means the residues read and reported by US-align from the supplied coordinates, not the full biological sequence or the span of residue numbers. Missing coordinates and residue-selection options can therefore affect the scale.
- **Short chains.** The standard protein distance-scale floor is used through 21 residues. Empty or missing reported chain lengths are rejected rather than falling back to the aligned-pair count.
- **Output parsing.** Missing or inconsistent distance tables preserve the match-symbol fallback, so such pairs are not distance-weighted. If full chain lengths cannot be read, the pair is reported as failed; the pipeline's existing missing-pair policy assigns TM-score 0 when building the tree. Check failure warnings before interpreting clusters.
- **Benchmark interpretation.** The MSA rule and batching changes should not be credited as a single performance intervention. Alignment-only benchmarks exclude downstream MSA and representative-selection costs; alignment quality and representative stability need separate comparisons.

### Running the tests

```bash
python -m pip install --no-deps --no-build-isolation -e ".[test]"
pytest
```

The reference Conda environment already includes the pinned test tools. Tests that need USalign, BeEM or TreeCluster.py are skipped when those tools are not on `PATH`; release validation should run with all three installed.

### References

BibTeX-formatted references for this project and the associated references can be found in [colabalign.bib](colabalign.bib) and [associated-references.bib](associated-references.bib).

- **ColabAlign (paper)**: <https://doi.org/10.1101/2025.10.06.677802>
- **ColabAlign (code only)**: <https://doi.org/10.5281/ZENODO.14169501>
- **BioPython**: <https://doi.org/10.1093/bioinformatics/btp163>
- **US-align**: <https://doi.org/10.1038/s41592-022-01585-1>
- **TreeCluster**: <https://doi.org/10.1371/journal.pone.0221068>
- **MUSTANG**: <https://doi.org/10.1002/prot.20921>
- **MView**: <https://doi.org/10.1093/bioinformatics/14.4.380>
- **BeEM**: <https://doi.org/10.1186/s12859-023-05388-9>

---

> [!NOTE]
> Permission to use, copy, modify, and distribute this program for any purpose, with or without fee, is hereby granted, provided that the notices on the head, the reference information, and this copyright notice appear in all copies or substantial portions of the Software. It is provided "as is" without express or implied warranty.
