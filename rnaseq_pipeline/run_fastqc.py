import os
import subprocess as sp

# FastQC
def run_fastqc(df, threads):
    """
    Run FastQC on a set of FASTQ files.

    Executes FastQC in quiet mode across the given FASTQ files, writing all
    reports to the "FastQC" directory (created if it doesn't already exist).

    Args:
        df (pandas.DataFrame): Sample sheet containing "name", "fastq1", and
            "fastq2" columns for each sample to process.
        threads (int): Number of threads to pass to FastQC for parallel
            processing of files.

    Raises:
        subprocess.CalledProcessError: If the fastqc command exits with a
            non-zero status.
    """
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)

    fastq_files = df[["fastq1", "fastq2"]].values.flatten().tolist()
    fastq_string = " ".join(str(f) for f in fastq_files)

    fastqc_command = (
        f"fastqc -q -t {threads} -o FastQC/ {fastq_string}"
    )

    sp.run(fastqc_command, shell=True, check=True)
