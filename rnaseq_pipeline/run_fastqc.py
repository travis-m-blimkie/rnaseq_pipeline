import os
import subprocess as sp

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
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)
    fastqc_command = (
        f"fastqc -q -t {threads} -o FastQC/ {all_fastq_files}"
    )
    sp.run(fastqc_command, shell=True, check=True)
