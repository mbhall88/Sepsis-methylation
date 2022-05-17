rule download_transcriptome:
    output:
        fasta=data_dir / "reference/GRCh38.rna.fa.gz",
    container:
        containers["base"]
    log:
        logs_dir / "download_transcriptome.log",
    resources:
        time="30m",
    params:
        url=transcriptome_url,
    shell:
        "wget -O {output.fasta} {params.url} 2> {log}"


rule index_transcriptome:
    input:
        fasta=rules.download_transcriptome.output.fasta,
    output:
        fasta=data_dir / "reference/transcriptome_reference.fa",
        fai=data_dir / "reference/transcriptome_reference.fa.fai",
    log:
        logs_dir / "index_transcriptome.log",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    container:
        containers["metacompore_python"]
    script:
        str(scripts_dir / "index_transcriptome.py")
