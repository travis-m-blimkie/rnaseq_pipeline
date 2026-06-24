# rnaseq_pipeline

Python implementation of a short read RNA-Seq pipeline, consisting of the following steps:

1. FastQC for overall read quality
2. Alignment to a reference genome with STAR
3. Gene-level read quantification with HTSeq
4. Results summary with MultiQC
5. Optional: Convert BAM files to CRAM format to reduce disk space usage

## Features

- Simple installation via Conda/Mamba and pip
- Support for multi-threading FastQC, STAR, and (optional) CRAM conversion
- Version information for each program in the pipeline is stored in a CSV file on completion

## Installation

```sh
# Clone the repository
git clone https://github.com/travis-m-blimkie/rnaseq_pipeline
cd rnaseq_pipeline

# Create the conda environment with all dependencies
conda env create -f environment.yml
```

## Example usage

Create a folder for your project, with a directory for the fastq files, genome index, and the sample table. The sample sheet should be a CSV file with the following content including exact column names:

| name | fastq1 | fastq2 |
|------|--------|--------|
| sample1 | Fastq/sample1_R1.fastq.gz | Fastq/sample1_R2.fastq.gz |
| sample2 | Fastq/sample2_R1.fastq.gz | Fastq/sample2_R2.fastq.gz |
| sample3 | Fastq/sample3_R1.fastq.gz | Fastq/sample3_R2.fastq.gz |

<br>

Then run the pipeline as follows:

```sh
conda activate rnaseq_pipeline

rnaseq_pipeline \
    --sample_sheet samples.csv \
    --genome_dir Genome \
    --threads 8 \
    --strand reverse \
    --gtf Genome/genome.gtf \
    --cram \
    --fasta Genome/genome.fa

conda deactivate
```

## Notes and limitations

- Only supports paired-end short read RNA-Seq data
- The genome index needs to be created separately with a compatible STAR version. We recommend installing `rnaseq_pipeline` first, then using the created conda environment to make the genome index
- HTSeq is hard-coded to run two processes simultaneously using GNU `parallel`
- Currently there is no support for passing additional command line arguments to individual programs; e.g. `htseq-count` is run with the default option for `--mode` (union)
- If performing BAM -> CRAM conversion with the `--cram` option and `--fasta` provided, the original BAM files are left in place and not deleted

## Links

- FastQC: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
- STAR: <https://alexdobin.github.io/STAR/>
- HTSeq: <https://htseq.readthedocs.io/en/latest/#>
- MultiQC: <https://seqera.io/multiqc/>
- Samtools: <https://www.htslib.org/>
- GNU parallel: <https://www.gnu.org/software/parallel/>

## Authors

Travis Blimkie

## AI disclosure

Claude and Posit Assistant was used to help write some sections of **rnaseq_pipeline**.
