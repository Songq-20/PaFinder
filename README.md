# PaFinder

PaFinder is a research tool for integrating candidate pathogen-related evidence in environmental bacterial MAGs. It is not a clinically validated pathogen classifier and does not produce pathogenicity probabilities.

## Features

- Pyrodigal ORF prediction
- VFDB Set A and SARG+ protein homology searches
- MGEdb nucleotide alignment
- VF/ARG–MGE genomic proximity and MGE loci
- skani ANI/AF reference matching
- Optional SVG/PNG genomic-neighborhood plots

## Requirements

Linux is recommended. The supported environment uses Python 3.12 and Conda. External programs are Pyrodigal, DIAMOND, NCBI BLAST+, and skani 0.4.0. Matplotlib is optional and is required only for `--plot-colocalization`.

## Installation

```bash
conda env create -f environment.yml
conda activate pafinder
python -m pip install --no-build-isolation -e .
# optional plotting
python -m pip install -e '.[plot]'
```

## Database preparation

The external database root passed with `--db` must contain normalized databases and a pathogen reference manifest:

```text
PaFinderDB/
└── build/databases/v1/
    ├── vfdb/{reference.faa,metadata.tsv}
    ├── sarg_plus/{reference.faa,metadata.tsv}
    ├── mgedb/{reference.fna,metadata.tsv}
    └── pathogen_references/{reference_manifest.tsv,skani_db}
```

Original databases are read-only. Build and validate pathogen references on the Linux/HPC host where FASTA files and skani are available.

## Quick start

```bash
pafinder bac -i MAGs/MAG_001.fna -o results --db /data/PaFinderDB -t 16
pafinder bac -i MAGs/ -o results --db /data/PaFinderDB -t 16
pafinder bac -i MAGs/ -o results --db /data/PaFinderDB -t 16 --skip-skani
pafinder bac -i MAGs/ -o results --db /data/PaFinderDB -t 16 --plot-colocalization
```

Use `--dry-run` to inspect the configured paths without running analysis.

## Output structure

```text
results/
├── summary.tsv
├── summary_stats.tsv
├── mags/<MAG_ID>/
│   ├── summary.tsv
│   ├── summary_stats.tsv
│   ├── evidence/{gene_summary.tsv,vf_hits.tsv,arg_hits.tsv,mge_hits.tsv,mge_loci.tsv,vf_mge_proximity.tsv,arg_mge_proximity.tsv,pathogen_matches.tsv}
│   ├── figures/
│   └── work/{phase3a/,phase3b/,skani_hits.raw.tsv}
└── qc/{mag_status.tsv,run_manifest.json}
```

The root summary has one row per input MAG. Evidence tables preserve raw alignments. Coordinates are 1-based inclusive. `mge_loci.tsv` contains overlapping query-coordinate loci and is not evidence of a complete mobile element.

## Evidence interpretation

VF and ARG homologs are sequence similarities requiring functional review. MGE HSPs are alignments to MGEdb references; MGE loci are coordinate clusters, not confirmed complete elements. Proximity is a genomic association and is not proof of horizontal gene transfer. ANI and alignment fractions describe relatedness to references and do not establish pathogenicity. Multiple references and ambiguities are retained.

## Reproducibility

Keep the database manifest, reference fingerprints, command parameters, `qc/run_manifest.json`, and software versions with each run. Reproduce on the same database snapshot and Conda environment.

## Testing

```bash
python -m compileall src
python -m pytest -q
```

Run a real MAG smoke test on Linux/HPC after installing Pyrodigal, DIAMOND, BLAST+, skani, and the normalized databases.

## Limitations

The first release targets bacterial MAGs. RNA viruses are not supported. Homology is not functional confirmation, database coverage is incomplete, and the software does not output a validated pathogenicity risk score.

## Citation / License

Citation and redistribution licensing information will be added when the project establishes the authoritative references and license terms. No DOI or license is asserted here.
