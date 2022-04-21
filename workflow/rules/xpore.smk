xpore_img = "docker://quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"


# we have to do an extra eventalign because xpore needs read_index instead of read_name
rule xpore_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.f5c_index.output.index,
        bam=rules.alignmemt_postfilter.output.bam,
        fasta=rules.index_transcriptome.output.fasta,
        kmer_model="resources/f5c/r9.4_70bps.u_to_t_rna.5mer.template.model",
    output:
        tsv=join("results", module_name, rule_name, "{cond}_{rep}_data.tsv"),
        summary=join("results", module_name, rule_name, "{cond}_{rep}_summary.tsv"),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        f5c_container
    shell:
        "f5c eventalign {params.opt} -t {threads} --kmer-model {input.kmer_model} -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary}  > {output.tsv} 2> {log}"



rule xpore_dataprep:
    input:
        eventalign=rules.xpore_eventalign.output.tsv,
    output:
        data=multiext(
            f"results/{module_name}/{rule_name}/{{cond}}_{{rep}}/data",
            ".index",
            ".json",
            ".log",
            ".readcount",
        ),
        idx=join("results", module_name, rule_name, "{cond}_{rep}", "eventalign.index"),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
        outdir=lambda wildcards, output: Path(output.idx).parent,
    resources:
        mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        xpore_img
    shell:
        """
        xpore dataprep {params.opt} --eventalign {input.eventalign} \
            --n_processes {threads} --out_dir {params.outdir} 2> {log}
        """



rule xpore_config:
    input:
        jsons=expand(
            f"results/{module_name}/xpore_dataprep/{{cond}}_{{rep}}/data.json",
            cond=condition_list,
            rep=replicates_list,
        ),
    output:
        configuration=join("results", module_name, rule_name, "config.yaml"),
    params:
        readcount_min=config[rule_name]["readcount_min"],
        readcount_max=config[rule_name]["readcount_max"],
        outdir=lambda wildcards, output: Path(output.configuration).parent.parent
        / "xpore_diffmod",
    log:
        join("logs", module_name, rule_name + ".log"),
    conda:
        f"../envs/xpore_config.yaml"
    script:
        f"../scripts/xpore_config.py"


rule xpore_diffmod:
    input:
        configuration=rules.xpore_config.output.configuration,
    output:
        table=join("results", module_name, rule_name, "diffmod.table"),
        log=join("results", module_name, rule_name, "diffmod.log"),
    log:
        join("logs", module_name, rule_name + ".log"),
    threads: get_threads(config, rule_name)
    resources:
        mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        xpore_img
    shell:
        "xpore diffmod --config {input.configuration} --n_processes {threads} 2> {log}"
