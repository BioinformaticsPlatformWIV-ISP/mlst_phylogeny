# MLST phylogeny
Script to construct phylogenies from tabular sequence typing output files.

----

## INSTALLATION

This script requires Python 3 and has only one external dependency: pandas, which can be installed with the 
following command:
```
pip install pandas
```

## USAGE

```
usage: mlst_phylo.py [-h] [--input-tsv INPUT_TSV] --output-matrix OUTPUT_MATRIX --output-dists OUTPUT_DISTS [--min-perc-loci MIN_PERC_LOCI] [--min-perc-samples MIN_PERC_SAMPLES]

options:
  -h, --help            show this help message and exit
  --input-tsv INPUT_TSV
  --output-matrix OUTPUT_MATRIX
                        Filtered allele matrix (TSV)
  --output-dists OUTPUT_DISTS
                        Pairwise distance matrix (TSV)
  --min-perc-loci MIN_PERC_LOCI
                        Minimum percentage of loci that should be present in a dataset
  --min-perc-samples MIN_PERC_SAMPLES
                        Minimum percentage of datasets where loci should be present
```

### Example command
 
```
mlst_phylo.py \
  --input-tsv typing-cgmlst-S1.tsv \
  --input-tsv typing-cgmlst-S2.tsv \
  --input-tsv typing-cgmlst-S3.tsv \
  --input-tsv typing-cgmlst-S4.tsv \
  --output-matrix matrix.tsv \
  --output-dists dists.tsv
```

## INPUT FORMAT

The tool requires several (at least three) TSV files with the following info:

| Locus      | Allele | % Identity | HSP/Locus length | Type |
|------------|--------|------------|------------------|------|
| SC0831     | 1      | 100.00     | 129/129          | DNA  |
| SEN0401    | 10     | 100.00     | 117/117          | DNA  |
| SPAB_04503 | 4      | 100.00     | 108/108          | DNA  |
| STM0834    | -      | -          | -                | DNA  |


The `Type` column is optional and is not used by the script.

Missing loci should be marked with `-` in the `Allele` column.

Multi-hits should be marked with `?` in the `Allele` column.

## TESTDATA

There are test datasets provided in the `testdata` directory.
These were generated using the *Salmonella* cgMLST scheme.
