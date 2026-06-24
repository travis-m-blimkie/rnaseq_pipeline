import argparse
import os
import subprocess as sp
from glob import glob
from re import search

import pandas as pd

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
    help="Stranded option for HTSeq"
)

parser.add_argument(
    "--gtf",
    type=str,
    help="Path to GTF for HTSeq"
)

parser.add_argument(
    "--cram",
    action="store_true",
    default=True,
    help="Enable CRAM output (default: True)"
)

parser.add_argument(
    "--fasta",
    type=str,
    help="Path to fasta for Samtools CRAM conversion"
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

# FastQC
def run_fastqc(all_fastq_files, threads):
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)
    fastqc_command = (
        f"fastqc -q -t {threads} -o FastQC/ {all_fastq_files}"
    )
    sp.run(fastqc_command, shell=True, check=True)

# STAR
def run_star(df, genome_dir, threads):
    os.makedirs(os.path.dirname("STAR/"), exist_ok=True)

    for _, row in df.iterrows():
        name = row["name"]
        fastq1 = row["fastq1"]
        fastq2 = row["fastq2"]

        out_prefix = f"STAR/{name}_"

        star_cmd = (
            f"STAR "
            f"--runThreadN {threads} "
            f"--runMode alignReads "
            f"--limitBAMsortRAM 48000000000 "
            f"--genomeLoad LoadAndKeep "
            f"--genomeDir {genome_dir} "
            f"--readFilesIn {fastq1} {fastq2} "
            f"--readFilesCommand zcat "
            f"--outFileNamePrefix {out_prefix} "
            f"--outSAMtype BAM SortedByCoordinate"
        )

        sp.run(star_cmd, shell=True, check=True)

    bam_files = glob("STAR/*.bam")
    for bam_file in bam_files:
        os.rename(
            bam_file,
            bam_file.replace("_Aligned.sortedByCoord.out", ""),
        )

    sp.run(
        f"STAR --genomeDir {genome_dir} --genomeLoad Remove",
        shell=True,
        check=True,
    )

# HTSeq
def run_htseq(df, strand, gtf_file):
    os.makedirs("HTSeq", exist_ok=True)

    # Build one htseq-count command per sample
    commands = []
    for _, row in df.iterrows():
        name = row["name"]
        htseq_cmd = (
            f"htseq-count "
            f"-s {strand} "
            f"-a 10 "
            f"-f bam "
            f"-r pos "
            f"STAR/{name}.bam "
            f"{gtf_file} "
            f"> HTSeq/{name}.count"
        )
        commands.append(htseq_cmd)

    # Feed all commands to GNU parallel, running two at a time
    joined_cmds = "\n".join(commands)
    sp.run(
        ["parallel", "-j", "2"],
        input=joined_cmds,
        text=True,
        shell=False,
        check=True
    )

# MultiQC
def run_multiqc():
    multiqc_cmd = "multiqc -f -o MultiQC FastQC STAR HTSeq"
    sp.run(multiqc_cmd, shell=True, check=True)

# Samtools/CRAM
def run_samtools(df, threads, fasta):
    for _, row in df.iterrows():
        name = row["name"]

        cram_cmd = (
            f"samtools view "
            f"--threads {threads} "
            f"-C "
            f"-T {fasta} "
            f"-o STAR/{name}.cram "
            f"STAR/{name}.bam"
        )
        sp.run(cram_cmd, shell=True, check=True)

# Version information
def run_versions(run_cram):
    version_fastqc = (
        sp.check_output("fastqc --version", shell=True)
        .decode("utf-8")
        .replace("\n", "")
        .replace("FastQC v", "")
    )

    version_star = (
        sp.check_output("STAR --version", shell=True)
        .decode("utf-8")
        .replace("\n", "")
    )

    version_htseq_raw = sp.check_output(
        "htseq-count --help | grep version",
        shell=True,
        stderr=sp.STDOUT,
    ).decode("utf-8")
    version_htseq = search(
        "[0-9]\\.[0-9]{1,2}\\.[0-9]{1,2}",
        version_htseq_raw,
    ).group(0)

    version_multiqc_raw = (
        sp.check_output("multiqc --version", shell=True)
        .decode("utf-8")
        .replace("\n", "")
    )
    version_multiqc = search(
        "[0-9]{1,2}\\.[0-9]{1,2}",
        version_multiqc_raw,
    ).group(0)

    if not run_cram:
        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC"],
            "Version": [version_fastqc, version_star, version_htseq, version_multiqc],
        }
    else:
        version_samtools_raw = sp.check_output(
            "samtools --version | grep samtools",
            shell=True,
        ).decode("utf-8").replace("\n", "")
        version_samtools = search(
            "[0-9]{1,2}\\.[0-9]{1,2}",
            version_samtools_raw,
        ).group(0)

        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC", "Samtools"],
            "Version": [
                version_fastqc,
                version_star,
                version_htseq,
                version_multiqc,
                version_samtools,
            ],
        }

    program_version_df = pd.DataFrame(program_version_dict)
    program_version_df.to_csv("version_information.csv", index=False)

# Run the functions
# run_fastqc(fastq_string, threads)
# run_star(df, genome_dir, threads)
run_htseq(df, strand, gtf)
run_multiqc()
if cram:
    run_samtools(df, threads, fasta)
run_versions(cram)
