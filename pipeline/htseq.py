# HTSeq
def run_htseq(df, strand, gtf_file):
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
