import os
import subprocess as sp
from glob import glob

# STAR
def run_star(df, genome_dir, threads):
    """
    Align FASTQ files to a reference genome using STAR.

    Iterates over each row in the sample sheet, running STAR in
    alignReads mode on the corresponding fastq1/fastq2 files, with the
    genome kept loaded in shared memory ("LoadAndKeep") across samples for
    efficiency. Output BAM files are sorted by coordinate and written to the
    "STAR" directory, then renamed from STAR's default
    "<name>_Aligned.sortedByCoord.out.bam" naming to "<name>.bam". Once all
    samples have been processed, the shared genome index is removed from
    memory.

    Args:
        df (pandas.DataFrame): Sample sheet containing "name", "fastq1", and
            "fastq2" columns for each sample to align.
        genome_dir (str): Path to the STAR genome index directory.
        threads (int): Number of threads to pass to STAR for alignment.

    Raises:
        subprocess.CalledProcessError: If any STAR command (alignment or
            genome unload) exits with a non-zero status.
    """
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

        sp.run(star_cmd, shell=True, check=True)

    bam_files = glob("STAR/*.bam")
    for bam_file in bam_files:
        os.rename(
            bam_file,
            bam_file.replace("_Aligned.sortedByCoord.out", ""),
        )

    sp.run(
        f"STAR --genomeDir {genome_dir} --genomeLoad Remove",
        shell=True,
        check=True,
    )

    files_to_remove = [
        "Aligned.out.sam",
        "Log.final.out",
        "Log.out",
        "Log.progress.out",
        "SJ.out.tab"
    ]

    for filename in files_to_remove:
        if os.path.exists(filename):
            os.remove(filename)
