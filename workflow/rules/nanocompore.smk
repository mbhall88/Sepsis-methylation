from pathlib import Path


rule nanocompore_eventalign_collapse:
    input:
        tsv=rules.f5c_eventalign.output.tsv,
    output:
        tsv=results / "nanocompore/{sample}/out_eventalign_collapse.tsv",
        idx=results / "nanocompore/{sample}/out_eventalign_collapse.tsv.idx",
    log:
        logs_dir / "nanocompore_eventalign_collapse/{sample}.log",
    threads: 4
    params:
        opt="",
        outdir=lambda wildcards, output: Path(output.tsv).parent,
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
        time="2d",
    container:
        containers["nanocompore"]
    shell:
        """
        nanocompore eventalign_collapse -t {threads} {params.opt} --overwrite \
            -i {input.tsv} -o {params.outdir} &> {log}
        """


rule nanocompore_sampcomp:
    input:
        control_tsv=results / f"nanocompore/{CTRL}/out_eventalign_collapse.tsv",
        test_tsv=rules.nanocompore_eventalign_collapse.output.tsv,
        fasta=rules.index_transcriptome.output.fasta,
    output:
        res_tsv=results / "nanocompore/{sample}/outnanocompore_results.tsv",
        shift_tsv=results / "nanocompore/{sample}/outnanocompore_shift_stats.tsv",
        res_db=results / "nanocompore/{sample}/outSampComp.db",
    wildcard_constraints:
        sample="|".join(TESTS),  # dont use control sample in {sample} wildcard
    log:
        logs_dir / "nanocompore_sampcomp/{sample}.log",
    threads: 4
    params:
        opt=" ".join(
            [
                "--max_invalid_kmers_freq 0.2 ",
                "--min_coverage 30 ",
                "--downsample_high_coverage 5000 ",
                "--min_ref_length 100 ",
                "--comparison_methods GMM,KS ",
                "--sequence_context 2 ",
                "--sequence_context_weights harmonic ",
                "--pvalue_thr 0.01 ",
                "--logit",
            ]
        ),
        outdir=lambda wildcards, output: Path(output.res_tsv).parent,
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
        time="4h",
    container:
        containers["nanocompore"]
    shell:
        """
        nanocompore sampcomp -t {threads} {params.opt} -w -1 {input.control_tsv} \
            -2 {input.test_tsv} --label1 control --label2 test -f {input.fasta} \
            -o {params.outdir} &> {log}
        """


rule nanocompore_postprocess:
    input:
        res_tsv=rules.nanocompore_sampcomp.output.res_tsv,
        fasta=rules.index_transcriptome.output.fasta,
    output:
        results / "nanocompore/{sample}/nanocompore_results_GMM_context_0.tsv",
        results / "nanocompore/{sample}/nanocompore_results_GMM_context_2.tsv",
        results / "nanocompore/{sample}/nanocompore_results_KS_dwell_context_0.tsv",
        results / "nanocompore/{sample}/nanocompore_results_KS_dwell_context_2.tsv",
        results / "nanocompore/{sample}/nanocompore_results_KS_intensity_context_0.tsv",
        results / "nanocompore/{sample}/nanocompore_results_KS_intensity_context_2.tsv",
    wildcard_constraints:
        sample="|".join(TESTS),  # dont use control sample in {sample} wildcard
    log:
        logs_dir / "nanocompore_postprocess/{sample}.log",
    params:
        opt={"p_val_lim": 0.01, "quantile_lim": 0.5, "min_distance": 9},
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
        time="20m",
    container:
        containers["metacompore_python"]
    script:
        str(scripts_dir / "nanocompore_postprocess.py")
