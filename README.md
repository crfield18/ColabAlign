# ColabAlign

## Fast pairwise protein secondary structure comparisons using multiprocessing

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/crfield18/ColabAlign/blob/main/colabalign.ipynb) [![DOI](https://zenodo.org/badge/788453062.svg)](https://doi.org/10.5281/zenodo.14169501)

This notebook performs pairwise protein structural alignments using the [_US-align_](https://zhanggroup.org/US-align/) algorithm by [_Zhang et al. (2022)_](https://doi.org/10.1038/s41592-022-01585-1) (an updated version of TM-align by [_Zhang and Skolnick (2005)_](https://doi.org/10.1093/nar/gki524)), then constructs a phylogenetic tree using the UPGMA algorithm to visualise similarities.

Designed to run directly in Google Colab for ease-of-use and to remove any local hardware requirements. This implementation also includes multiprocessing support for dramatically increased performance over the base US-align program.

---

### Installation for local usage

ColabAlign.py is designed to work for Google Colab and running on local machines. A [YAML file](colabalign.yml) is provided for easy installation of dependencies in a Conda environment.

**On Linux (distro-dependent) and x86 Macs (i.e. pre-M1)**, simply create an environment with:

`conda env create -f colabalign.yml`

**On ARM-based Macs (M1 and onwards)**, an extra flag is needed that allows x86-only scripts:

`conda env create --platform osx-64 -f colabalign.yml`

This also requires Rosetta 2 to be installed:

`softwareupdate --install-rosetta`

---

## References

BibTeX-formatted references for this project and the associated references can be found in [colabalign.bib](colabalign.bib) and [associated-references.bib](associated-references.bib).

- **ColabAlign**: <https://doi.org/10.5281/ZENODO.14169501>
- **BioPython**: <https://doi.org/10.1093/bioinformatics/btp163>
- **US-align**: <https://doi.org/10.1038/s41592-022-01585-1>
- **TreeCluster**: <https://doi.org/10.1371/journal.pone.0221068>
- **MUSTANG**: <https://doi.org/10.1002/prot.20921>
- **MView**: <https://doi.org/10.1093/bioinformatics/14.4.380>

---

> [!IMPORTANT]
> Permission to use, copy, modify, and distribute this program for any purpose, with or without fee, is hereby granted, provided that the notices on the head, the reference information, and this copyright notice appear in all copies or substantial portions of the Software. It is provided "as is" without express or implied warranty.
