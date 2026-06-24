import os
import subprocess as sp

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
