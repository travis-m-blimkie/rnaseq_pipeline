# STAR
def run_star(df, genome_dir, threads):
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
