import subprocess as sp

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
