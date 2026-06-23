import os
import subprocess as sp
import pandas as pd
from glob import glob
from re import search
import argparse

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
    "--do_cram",
    action="store_true",
    default=True,
    help="Enable CRAM output (default: True)"
)

args = parser.parse_args()

## Define input parameters
sample_sheet = args.sample_sheet
genome_dir = args.genome_dir
threads = args.threads
strand = args.strand
gtf = args.gtf
do_cram = args.do_cram


df = pd.read_csv(sample_sheet)

fastq_files = df[["fastq1", "fastq2"]].values.flatten().tolist()
fastq_string = " ".join(str(f) for f in fastq_files)


## FastQC
def run_fastqc(all_fastq_files, threads):
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)
    fastqc_command = "fastqc -q -t " + str(threads) + " -o FastQC/ " + all_fastq_files
    sp.run(fastqc_command, shell=True)


## STAR
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

        sp.run(star_cmd, shell=True)

    bam_files = glob("STAR/*.bam")
    for b in bam_files:
        os.rename(b, b.replace("_Aligned.sortedByCoord.out", ""))

    sp.run(f"STAR --genomeDir {genome_dir} --genomeLoad Remove", shell=True)
    sp.run("rmdir _STARtmp/", shell=True)
    sp.run("rm Aligned.out.sam Log.out Log.final.out Log.progress.out SJ.out.tab", shell=True)



## HTSeq
def run_htseq(df, strand, gtf_file):
    os.makedirs(os.path.dirname("HTSeq/"), exist_ok=True)

    for _, row in df.iterrows():
        name = row["name"]

        htseq_cmd = (
            f"htseq-count "
            f"-s {strand} "
            f"-a 10 "
            f"-f bam "
            f"-r pos "
            f"-c HTSeq/{name}.count "
            f"STAR/{name}.bam "
            f"{gtf_file} "
        )
        print(htseq_cmd)


## MultiQC
def run_multiqc():
    multiqc_cmd = "multiqc -f -o MultiQC FastQC STAR HTSeq"
    sp.run(multiqc_cmd, shell=True)


## Samtools/CRAM
# def run_samtools(fasta_file):
#     bam_files = glob("STAR/*.bam")
#     for b in bam_files:
#         cram_name = b.replace("bam", "cram")
#         cram_command = "samtools view --threads " + threads + " -C " + b + " -T " + fasta_file + " > " + cram_name

#         if exists(cram_name):
#             print("Skipping, CRAM file already exists for", cram_name)
#         else:
#             sp.run(cram_command, shell=True)


## Version information
def run_versions(run_cram):
    version_fastqc = sp.check_output("fastqc --version", shell=True).decode("utf-8").replace("\n", "").replace("FastQC v", "")

    version_star = sp.check_output("STAR --version", shell=True).decode("utf-8").replace("\n", "")

    version_htseq_raw = sp.check_output("htseq-count --help | grep version", shell=True, stderr=sp.STDOUT).decode("utf-8")
    version_htseq = search("[0-9]\\.[0-9]{1,2}\\.[0-9]{1,2}", version_htseq_raw).group(0)

    version_multiqc_raw = sp.check_output("multiqc --version", shell=True).decode("utf-8").replace("\n", "")
    version_multiqc = search("[0-9]{1,2}\\.[0-9]{1,2}", version_multiqc_raw).group(0)

    if not(run_cram):
        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC"],
            "Version": [version_fastqc, version_star, version_htseq, version_multiqc]
        }
    elif run_cram :
        version_samtools_raw = sp.check_output("samtools --version | grep samtools", shell=True).decode("utf-8").replace("\n", "")
        version_samtools = search("[0-9]{1,2}\\.[0-9]{1,2}", version_samtools_raw).group(0)

        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC", "Samtools"],
            "Version": [version_fastqc, version_star, version_htseq, version_multiqc, version_samtools]
        }
    
    program_version_df = pd.DataFrame(program_version_dict)
    program_version_df.to_csv("version_information.csv", index=False)


## Run all the functions
run_fastqc(fastq_string, threads)
#run_star(input_dir=my_fastq_dir, genome_dir=my_genome_dir)
#run_htseq(gtf_file=my_gtf_file)
#run_multiqc()
#if do_cram:
#    run_samtools(fasta_file=my_fasta_file)
run_versions(run_cram=do_cram)
