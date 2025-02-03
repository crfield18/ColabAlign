# ColabAlign

## Fast pairwise protein secondary structure comparisons using multiprocessing

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/crfield18/ColabAlign/blob/main/colabalign.ipynb) [![DOI](https://zenodo.org/badge/788453062.svg)](https://doi.org/10.5281/zenodo.14169501)

This notebook performs pairwise protein structural alignments using the [_US-align_](https://zhanggroup.org/US-align/) algorithm by [_Zhang et al. (2022)_](https://doi.org/10.1038/s41592-022-01585-1) (an updated version of TM-align by [_Zhang and Skolnick (2005)_](https://doi.org/10.1093/nar/gki524)), then constructs a phylogenetic tree using the UPGMA algorithm to visualise similarities.

Designed to run directly in Google Colab for ease-of-use and to remove any local hardware requirements. This implementation also includes multiprocessing support for dramatically increased performance over the base US-align program.

Made possible with [US-align](https://github.com/pylelab/USalign) and [Biopython](https://biopython.org). Based on [_mTM-align_](http://yanglab.nankai.edu.cn/mTM-align/) by [_Dong et al. (2018)_](https://doi.org/10.1093/nar/gky430).
