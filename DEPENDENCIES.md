# ColabAlign v1.2.0 dependency baseline

The v1.2.0 package uses exact direct dependency versions rather than resolving whatever versions happen to be current at installation time. The complete Conda locks also freeze transitive packages, binary builds and archive hashes. These pins describe the new release, not the environment used for the manuscript's original candidate selection.

## Reference versions

| Component | Version |
|---|---|
| Python | 3.12.16 |
| NumPy | 2.5.4 |
| pandas | 3.1.1 |
| Biopython | 1.88 |
| SciPy | 1.18.1 |
| tqdm | 4.70.1 |
| TreeCluster | 1.0.5 |
| US-align | 20241201 |
| BeEM (`beem-bio`) | 1.0.1 |
| pip | 26.2.2 |
| setuptools | 84.2.1 |
| wheel | 0.48.0 |
| pytest | 9.2.0 |
| build (`python-build` in Conda) | 1.6.1 |
| Konda (Colab bootstrap only) | 0.1.0 |

`pyproject.toml` pins the direct Python runtime, build-system, test and release dependencies. `colabalign.yml` and the Conda recipe pin the corresponding versions plus the interpreter and external tools. Package metadata supports the Python 3.12 series; the reference environment and Conda recipe specify the exact interpreter patch version. Other Python minor versions are not claimed to be tested by this release.

Biopython 1.88 includes the UPGMA clade-height correction introduced in 1.85. Reproducing a historical dendrogram generated with older Biopython may require that historical environment rather than this release baseline. See the [Biopython correction](https://github.com/biopython/biopython/pull/3951) and the manuscript's provenance discussion.

## Exact Conda environments

The explicit locks in `locks/` include every resolved Conda dependency, its build, download URL and archive MD5. They are generated using micromamba 2.10.0 with strict channel priority, `conda-forge` before `bioconda`, and no defaults channel.

- `conda-linux-64.lock` is the installed Linux x86-64 reference environment.
- `conda-osx-64.lock` is a solver-validated Intel macOS environment.
- `conda-osx-arm64.lock` is a solver-validated native Apple Silicon environment.

The macOS locks were solved with `CONDA_OVERRIDE_OSX=11.0` to model macOS 11 or later. Solving these environments on Linux is not a substitute for executing the test suite on each Mac architecture.

Create the environment with the appropriate lock, then install ColabAlign from the exact release tag or commit being tested:

```bash
conda create -n colabalign --file locks/conda-linux-64.lock
conda activate colabalign
python -m pip install --no-deps --no-build-isolation .
python -m pip check
python -m pytest -q
python -m build --no-isolation
```

`micromamba create -n colabalign --file ...` can be used instead of `conda create`. The locks are platform-specific; do not use a Linux lock on a Mac. They intentionally include the build and test tools so the same environment can validate and build the release.

`conda env create -f colabalign.yml` keeps direct versions fixed but resolves transitive versions and build variants afresh. Use the explicit lock for exact reproduction, not just the YAML.

## Exact Python-only environment

`locks/requirements-py312.txt` freezes the Python runtime, test, release and transitive dependencies with SHA-256 archive hashes. It does not install Python itself or the US-align and BeEM executables.

```bash
# Start with a clean Python 3.12.16 virtual environment.
python -m pip install --require-hashes -r locks/requirements-py312.txt
python -m pip install --no-deps --no-build-isolation .
```

Install US-align 20241201 and BeEM 1.0.1 separately, or use one of the full Conda locks. Do not apply the pip lock over a working Conda lock unless deliberately replacing its package provenance.

The public Colab notebook follows `main`; `benchmarking` is for local testing only. The notebook pins its Konda bootstrap to 0.1.0 and creates a separate `colabalign` analysis environment from the explicit Linux lock instead of updating Google's base environment. Installation and analysis both use that environment. This freezes analysis dependencies but not Google's notebook kernel, hosted infrastructure or the Miniconda installer fetched by Konda. The hosted notebook has not been exercised by the local Linux package tests.

## Updating the pins

Update `pyproject.toml`, `colabalign.yml` and `conda-recipe/colabalign/meta.yaml` together. Re-solve the YAML using strict channel priority for every supported platform, export exact URLs and hashes to the explicit locks, and regenerate the Python lock:

```bash
uv pip compile pyproject.toml --extra test --extra release \
    --python-version 3.12 --universal --generate-hashes \
    --output-file locks/requirements-py312.txt
```

For macOS cross-solves, set `CONDA_OVERRIDE_OSX=11.0` and select the target platform. Do not invent build strings or hashes, or copy them between architectures. Recreate at least the native test environment from the resulting lock and run the complete suite, including tests requiring US-align, BeEM and TreeCluster.

`tests/test_dependency_pins.py` checks direct-pin agreement and that every direct Conda dependency appears in all three explicit locks with a matching version and a hash.

## Release checks still required

- Run the full test suite on native Intel and Apple Silicon macOS.
- Replace the existing all-zero source checksum in the Conda recipe with the actual SHA-256 of the `v1.2.0` tag archive once that tag exists. A source checksum cannot be derived from an unpublished tag.
- Build and validate the final Conda package from that archive before submitting it to Bioconda.
- Archive the released source, the lockfiles and benchmark provenance alongside the release.

The MSA distance-scale revision and v2 cache namespace are documented in [README.md](README.md#structure-informed-msa-distance-weights). Dependency pinning does not establish equivalence with historical algorithms or environments.
