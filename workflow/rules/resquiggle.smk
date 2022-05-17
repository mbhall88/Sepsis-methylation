
rule f5c_index:
    input:
        fast5_dir=get_fast5_dir,
        fastq=rules.merge_fastq.output.fastq,
    output:
        indices=multiext(
            rules.merge_fastq.output.fastq,
            ".index",
            ".index.fai",
            ".index.gzi",
            ".index.readdb",
        ),
    log:
        logs_dir / "f5c_index/{sample}.log",
    threads: 4
    resources:
        partition="gpgpu",
        slurm="gres='gpu:1' qos=gpgpumdhs",
    params:
        opt="--iop 4",
    container:
        containers["f5c"]
    shell:
        "f5c index {params.opt} -t {threads} -d {input.fast5_dir} {input.fastq} 2> {log}"


rule f5c_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.f5c_index.output.indices,
        bam=rules.alignmemt_postfilter.output.bam,
        fasta=rules.index_transcriptome.output.fasta,
        kmer_model="resources/f5c/r9.4_70bps.u_to_t_rna.5mer.template.model",
    output:
        tsv=results / "eventalign/{sample}_data.tsv",
        summary=results / "eventalign/{sample}_summary.tsv",
    log:
        logs_dir / "f5c_eventalign/{sample}.log",
    threads: 8
    params:
        opt="-x desktop-high --rna --samples --signal-index --print-read-names --scale-events --verbose 2",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * GB,
        partition="gpgpu",
        slurm="gres='gpu:1' qos=gpgpumdhs",
    container:
        containers["f5c"]
    shell:
        """
        f5c eventalign {params.opt} -t {threads} --kmer-model {input.kmer_model} \
            -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary} \
            > {output.tsv} 2> {log}
        """
