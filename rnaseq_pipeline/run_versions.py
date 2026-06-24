import subprocess as sp
from re import search
import pandas as pd

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

    version_parallel_raw = sp.run(['parallel', '--version'], capture_output=True, text=True)
    version_parallel_match = search(r'GNU parallel (\d{8})', version_parallel_raw.stdout)
    version_parallel = version_parallel_match.group(1)

    if not run_cram:
        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "parallel", "MultiQC"],
            "Version": [version_fastqc, version_star, version_htseq, version_parallel, version_multiqc],
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
            "Program": ["FastQC", "STAR", "HTSeq", "parallel", "MultiQC", "Samtools"],
            "Version": [
                version_fastqc,
                version_star,
                version_htseq,
                version_parallel,
                version_multiqc,
                version_samtools,
            ],
        }

    program_version_df = pd.DataFrame(program_version_dict)
    program_version_df.to_csv("version_information.csv", index=False)
