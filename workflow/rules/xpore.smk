from pathlib import Path

xpore_dir = results / "xpore"


# we have to do an extra eventalign because xpore needs read_index instead of read_name
rule xpore_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.f5c_index.output.indices,
        bam=rules.alignmemt_postfilter.output.bam,
        fasta=rules.index_transcriptome.output.fasta,
        kmer_model="resources/f5c/r9.4_70bps.u_to_t_rna.5mer.template.model",
    output:
        tsv=xpore_dir / "eventalign/{sample}_data.tsv",
        summary=xpore_dir / "eventalign/{sample}_summary.tsv",
    log:
        logs_dir / "xpore_eventalign/{sample}.log",
    threads: 8
    params:
        opt="-x desktop-high --rna --signal-index --scale-events --verbose 2",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * GB,
        partition="gpgpu",
        slurm="gres=gpu:1 qos=gpgpumdhs",
    container:
        containers["f5c"]
    shell:
        """
        f5c eventalign {params.opt} -t {threads} --kmer-model {input.kmer_model} \
            -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary} \
            > {output.tsv} 2> {log}
        """


rule xpore_dataprep:
    input:
        eventalign=rules.xpore_eventalign.output.tsv,
    output:
        data=multiext(
            str(xpore_dir / "dataprep/{sample}/data"),
            ".index",
            ".json",
            ".log",
            ".readcount",
        ),
        idx=xpore_dir / "{sample}/eventalign.index",
    log:
        logs_dir / "xpore_dataprep/{sample}.log",
    threads: 8
    params:
        opt="--readcount_max 1000000",
        outdir=lambda wildcards, output: Path(output.idx).parent,
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * GB,
        time="3h",
    container:
        containers["xpore"]
    shell:
        """
        xpore dataprep {params.opt} --eventalign {input.eventalign} \
            --n_processes {threads} --out_dir {params.outdir} 2> {log}
        """


diffmod_dir = xpore_dir / "diffmod"


rule xpore_config:
    input:
        control_json=xpore_dir / f"dataprep/{CTRL}/data.json",
        test_json=xpore_dir / "dataprep/{sample}/data.json",
    output:
        configuration=xpore_dir / "{sample}.config.yaml",
    wildcard_constraints:
        sample="|".join(TESTS),  # dont use control sample in {sample} wildcard
    params:
        readcount_min=15,
        readcount_max=1_000,
        outdir=lambda wildcards: diffmod_dir / wildcards.sample,
    resources:
        time="5m",
    log:
        logs_dir / "xpore_config/{sample}.log",
    conda:
        str(envs_dir / "xpore_config.yaml")
    script:
        str(scripts_dir / "xpore_config.py")


rule xpore_diffmod:
    input:
        configuration=rules.xpore_config.output.configuration,
    output:
        table=diffmod_dir / "{sample}/diffmod.table",
        log=diffmod_dir / "{sample}/diffmod.log",
    wildcard_constraints:
        sample="|".join(TESTS),  # dont use control sample in {sample} wildcard
    log:
        logs_dir / "xpore_diffmod/{sample}.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * GB,
        time="1h",
    container:
        containers["xpore"]
    shell:
        "xpore diffmod --config {input.configuration} --n_processes {threads} 2> {log}"
