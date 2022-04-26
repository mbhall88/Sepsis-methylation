
rule nanocompore_eventalign_collapse:
    input:
        tsv=rules.f5c_eventalign.output.tsv,
    output:
        outdir=directory(join("results", module_name, rule_name, "{cond}_{rep}")),
        tsv=join(
            "results",
            module_name,
            rule_name,
            "{cond}_{rep}",
            "out_eventalign_collapse.tsv",
        ),
        idx=join(
            "results",
            module_name,
            rule_name,
            "{cond}_{rep}",
            "out_eventalign_collapse.tsv.idx",
        ),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    # resources:
    #     mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        containers["nanocompore"]
    shell:
        "nanocompore eventalign_collapse -t {threads} {params.opt} --overwrite -i {input.tsv} -o {output.outdir} &> {log}"


rule nanocompore_sampcomp:
    input:
        control_tsv=expand(
            join(
                "results",
                module_name,
                "nanocompore_eventalign_collapse",
                "control_{rep}",
                "out_eventalign_collapse.tsv",
            ),
            rep=replicates_list,
        ),
        test_tsv=expand(
            join(
                "results",
                module_name,
                "nanocompore_eventalign_collapse",
                "test_{rep}",
                "out_eventalign_collapse.tsv",
            ),
            rep=replicates_list,
        ),
        fasta=rules.index_transcriptome.output.fasta,
    output:
        res_tsv=join("results", module_name, rule_name, "outnanocompore_results.tsv"),
        shift_tsv=join(
            "results", module_name, rule_name, "outnanocompore_shift_stats.tsv"
        ),
        res_db=join("results", module_name, rule_name, "outSampComp.db"),
    log:
        join("logs", module_name, rule_name, "sampcomp.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    # resources:
    #     mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        containers["nanocompore"]
    script:
        str(scripts_dir / "nanocompore_sampcomp.py")


rule nanocompore_postprocess:
    input:
        res_tsv=rules.nanocompore_sampcomp.output.res_tsv,
        fasta=rules.index_transcriptome.output.fasta,
    output:
        results / "final/{sample}/nanocompore_results_GMM_context_0.tsv",
        results / "final/{sample}/nanocompore_results_GMM_context_2.tsv",
        results / "final/{sample}/nanocompore_results_KS_dwell_context_0.tsv",
        results / "final/{sample}/nanocompore_results_KS_dwell_context_2.tsv",
        results / "final/{sample}/nanocompore_results_KS_intensity_context_0.tsv",
        results / "final/{sample}/nanocompore_results_KS_intensity_context_2.tsv",
    log:
        logs_dir / "nanocompore_postprocess/{sample}.log",
    params:
        opt={"p_val_lim": 0.01, "quantile_lim": 0.5, "min_distance": 9},
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    container:
        containers["metacompore_python"]
    script:
        str(scripts_dir / "nanocompore_postprocess.py")
