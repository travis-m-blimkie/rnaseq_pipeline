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
    """
    Run FastQC on a set of FASTQ files.

    Executes FastQC in quiet mode across the given FASTQ files, writing all
    reports to the "FastQC" directory (created if it doesn't already exist).

    Args:
        all_fastq_files (str): Space-separated string of paths to FASTQ
            files to run FastQC on.
        threads (int): Number of threads to pass to FastQC for parallel
            processing of files.

    Raises:
        subprocess.CalledProcessError: If the fastqc command exits with a
            non-zero status.
    """
    os.makedirs("FastQC", exist_ok=True)
    fastqc_command = (
        f"fastqc -q -t {threads} -o FastQC/ {all_fastq_files}"
    )
    sp.run(fastqc_command, shell=True, check=True)

# STAR
def run_star(df, genome_dir, threads):
    """
    Align FASTQ files to a reference genome using STAR.

    Iterates over each row in the sample sheet, running STAR in
    alignReads mode on the corresponding fastq1/fastq2 files, with the
    genome kept loaded in shared memory ("LoadAndKeep") across samples for
    efficiency. Output BAM files are sorted by coordinate and written to the
    "STAR" directory, then renamed from STAR's default
    "<name>_Aligned.sortedByCoord.out.bam" naming to "<name>.bam". Once all
    samples have been processed, the shared genome index is removed from
    memory.

    Args:
        df (pandas.DataFrame): Sample sheet containing "name", "fastq1", and
            "fastq2" columns for each sample to align.
        genome_dir (str): Path to the STAR genome index directory.
        threads (int): Number of threads to pass to STAR for alignment.

    Raises:
        subprocess.CalledProcessError: If any STAR command (alignment or
            genome unload) exits with a non-zero status.
    """
    os.makedirs("STAR", exist_ok=True)

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
    """
    Quantify gene counts from BAM files using htseq-count, run via GNU parallel.

    Builds one htseq-count command per sample in the sample sheet, reading
    each sample's BAM file ("STAR/<name>.bam") and writing counts to
    "HTSeq/<name>.count". All commands are passed to GNU parallel, which
    runs two samples concurrently.

    Args:
        df (pandas.DataFrame): Sample sheet containing a "name" column used
            to locate each sample's BAM file.
        strand (str): Strandedness setting passed to htseq-count's -s flag
            (e.g. "yes", "no", "reverse").
        gtf_file (str): Path to the GTF annotation file used for counting.

    Raises:
        subprocess.CalledProcessError: If GNU parallel exits with a
            non-zero status, e.g. because one or more htseq-count jobs
            failed.
    """
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
    """
    Run MultiQC to aggregate QC reports from FastQC, STAR, and HTSeq outputs.

    Searches the "FastQC", "STAR", and "HTSeq" directories for compatible
    log/report files and compiles them into a single summary report written
    to the "MultiQC" directory. Existing reports in the output directory are
    overwritten.

    Raises:
        subprocess.CalledProcessError: If the multiqc command exits with a
            non-zero status.
    """
    multiqc_cmd = "multiqc -f -o MultiQC FastQC STAR HTSeq"
    sp.run(multiqc_cmd, shell=True, check=True)

# Samtools/CRAM
def run_samtools(df, threads, fasta):
    """
    Convert per-sample BAM files to CRAM format using samtools.

    Iterates over each row in the sample sheet, converting the corresponding
    STAR-aligned BAM file ("STAR/<name>.bam") to a reference-compressed CRAM
    file ("STAR/<name>.cram") in the same directory.

    Args:
        df (pandas.DataFrame): Sample sheet containing a "name" column used
            to locate each sample's BAM file.
        threads (int): Number of threads to pass to samtools for compression.
        fasta (str): Path to the reference FASTA file used for CRAM
            reference-based compression (passed to samtools -T).

    Raises:
        subprocess.CalledProcessError: If any samtools command exits with a
            non-zero status.
    """
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
    """
    Record the installed versions of all pipeline tools to a CSV file.

    Queries the installed versions of FastQC, STAR, HTSeq, and MultiQC by
    calling each tool's version command and parsing the output. If
    `run_cram` is True, samtools' version is also queried and included,
    since CRAM conversion depends on it. The collected program/version
    pairs are written to "version_information.csv" in the current working
    directory.

    Args:
        run_cram (bool): Whether the pipeline includes a CRAM conversion
            step. If True, samtools is included in the version report.

    Raises:
        subprocess.CalledProcessError: If any version-check command exits
            with a non-zero status.
        AttributeError: If a version string cannot be parsed from a tool's
            output (i.e. the regex search returns no match).
    """
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
# run_htseq(df, strand, gtf)
run_multiqc()
if cram:
    run_samtools(df, threads, fasta)
run_versions(cram)
