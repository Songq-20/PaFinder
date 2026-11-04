# PaFinder v1

Phase 1 database standardization. Conda environment:

```bash
conda env create -f environment.yml
conda activate pafinder
python scripts/prepare_databases.py
```

Original databases are read-only inputs; normalized outputs are written to `build/databases/v1/`.

## Installation

```bash
conda env create -f environment.yml
conda activate pafinder
pip install .
```

Optional plotting dependency: `pip install '.[plot]'`.

## CLI status

The packaged `pafinder bac` entry point currently assembles existing Phase 3A/3B outputs into final TSV evidence files. It does not yet execute Pyrodigal, DIAMOND, BLASTN, or skani from a fresh MAG directory; those tools are available through the phase-specific scripts. Use `--dry-run` to inspect paths.
