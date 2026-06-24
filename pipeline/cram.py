# Samtools/CRAM
def run_samtools(df, threads, fasta):
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
