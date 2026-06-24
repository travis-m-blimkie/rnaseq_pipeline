import os
import subprocess as sp

# FastQC
def run_fastqc(all_fastq_files, threads):
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)
    fastqc_command = (
        f"fastqc -q -t {threads} -o FastQC/ {all_fastq_files}"
    )
    sp.run(fastqc_command, shell=True, check=True)
