import subprocess as sp

# MultiQC
def run_multiqc():
    multiqc_cmd = "multiqc -f -o MultiQC FastQC STAR HTSeq"
    sp.run(multiqc_cmd, shell=True, check=True)
