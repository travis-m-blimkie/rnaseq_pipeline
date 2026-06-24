import argparse
import os
import subprocess as sp
from glob import glob
from re import search

import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Pipeline parameters")

    parser.add_argument(
        "--sample_sheet",
        type=str,
        required=True,
        help="Path to the sample sheet CSV file"
    )

    parser.add_argument(
        "--genome_dir",
        type=str,
        required=True,
        help="Path to the genome directory"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use (default: 8)"
    )

    parser.add_argument(
        "--strand",
        type=str,
        default="reverse",
        help="Stranded option for HTSeq: yes, no, reverse (default)"
    )

    parser.add_argument(
        "--gtf",
        type=str,
        help="Path to the GTF file for HTSeq"
    )

    parser.add_argument(
        "--cram",
        action="store_true",
        help="Enable CRAM output"
    )

    parser.add_argument(
        "--fasta",
        type=str,
        help="Path to fasta for Samtools CRAM conversion. Required if --cram is specified"
    )

    args = parser.parse_args()

    # Define input parameters
    sample_sheet = args.sample_sheet
    genome_dir = args.genome_dir
    threads = args.threads
    strand = args.strand
    gtf = args.gtf
    cram = args.cram
    fasta = args.fasta

    df = pd.read_csv(sample_sheet)

    fastq_files = df[["fastq1", "fastq2"]].values.flatten().tolist()
    fastq_string = " ".join(str(f) for f in fastq_files)

    # Run the functions
    run_fastqc(fastq_string, threads)
    run_star(df, genome_dir, threads)
    run_htseq(df, strand, gtf)
    run_multiqc()
    if cram:
        run_samtools(df, threads, fasta)
    run_versions(cram)

if __name__ == "__main__":
    main()
