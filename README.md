# rnaseq_pipeline

Python implementation of a short read RNA-Seq pipeline, consisting of the following steps:

1. FastQC for overall read quality
2. Alignment to a specified reference genome with STAR
3. Gene-level read quantification with HTSeq
4. Results summary with MultiQC
5. Optional: Convert BAM files to CRAM format to reduce disk space usage

When completed, a CSV file is created to store version information for each program used by the pipeline.

## Additional requirements

- Assumes each of the programs is already installed and available (Python, FastQC, STAR, HTSeq, Samtools, MultiQC)
- Genome index should already be created with a compatible STAR version
- Only supports paired-end short read RNA-Seq data
- Currently no support for passing additional command line arguments to programs; e.g. `htseq-count` is run with the default option for `--mode` (union)
