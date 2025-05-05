# MLST phylogeny
Script to construct phylogenies from tabular sequence typing output files.

----

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


## TESTDATA

There are test datasets provided in the `testdata` directory.
These were generated using the *Salmonella* cgMLST scheme.
