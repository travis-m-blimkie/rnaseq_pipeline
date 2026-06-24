import subprocess as sp

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
