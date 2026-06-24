import subprocess as sp
from re import search
import pandas as pd

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
