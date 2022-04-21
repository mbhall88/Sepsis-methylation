
rule index_transcriptome:
    input:
        fasta=transcriptome_ref,
    output:
        fasta=join("results", module_name, rule_name, "transcriptome_reference.fa"),
        fai=join("results", module_name, rule_name, "transcriptome_reference.fa.fai"),
    log:
        join("logs", module_name, rule_name, "transcriptome_reference.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    container:
        containers["metacompore_python"]
    script:
        str(scripts_dir / "index_transcriptome.py")


rule minimap2_index:
    input:
        fasta=rules.index_transcriptome.output.fasta,
    output:
        idx=join("results", module_name, rule_name, "transcriptome_reference.mmi"),
    log:
        join("logs", module_name, rule_name, "name.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    conda:
        str(envs_dir / "aln_tools.yaml")
    shell:
        "minimap2 -t {threads} {params.opt} -d {output.idx} {input.fasta} 2> {log}"


rule minimap2_align:
    input:
        idx=rules.minimap2_index.output.idx,
        fastq=rules.merge_fastq.output.fastq,
    output:
        bam=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam")),
        bam_index=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam.bai")),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    conda:
        str(envs_dir / "aln_tools.yaml")
    shell:
        """
        minimap2 -t {threads} {params.opt} {input.idx} {input.fastq} 2> {log} | \
        samtools view -bh 2>> {log} | \
        samtools sort > {output.bam} 2>> {log} \
        && samtools index {output.bam} 2>> {log}
        """


rule alignmemt_prefilter:
    input:
        bam=rules.minimap2_align.output.bam,
        bam_index=rules.minimap2_align.output.bam_index,
    output:
        bam=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam")),
        bam_index=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam.bai")),
        reads_index=temp(
            join("results", module_name, rule_name, "{cond}_{rep}.bam.idx.gz")
        ),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    container:
        containers["pybiotools"]
    shell:
        "pyBioTools Alignment Filter {params.opt} -i {input.bam} -o {output.bam} --verbose &> {log}"


rule min_ref_coverage:
    input:
        reads_index_list=expand(
            join(
                "results",
                module_name,
                "alignmemt_prefilter",
                "{cond}_{rep}.bam.idx.gz",
            ),
            cond=condition_list,
            rep=replicates_list,
        ),
    output:
        ref_list=join("results", module_name, rule_name, "valid_references_list.txt"),
    log:
        join("logs", module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    container:
        containers["metacompore_python"]
    script:
        "../scripts/min_ref_coverage.py"


rule alignmemt_postfilter:
    input:
        bam=rules.alignmemt_prefilter.output.bam,
        bam_index=rules.alignmemt_prefilter.output.bam_index,
        ref_list=rules.min_ref_coverage.output.ref_list,
    output:
        bam=join("results", module_name, rule_name, "{cond}_{rep}.bam"),
        bam_index=join("results", module_name, rule_name, "{cond}_{rep}.bam.bai"),
        reads_index=join("results", module_name, rule_name, "{cond}_{rep}.bam.idx.gz"),
        selected_reads_fn=join(
            "results", module_name, rule_name, "{cond}_{rep}_selected_read_ids.txt"
        ),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    container:
        containers["pybiotools"]
    shell:
        """
        pyBioTools Alignment Filter {params.opt} -i {input.bam} \
            --select_ref_fn {input.ref_list} -o {output.bam} \
            -l {output.selected_reads_fn} --verbose &> {log}
        """
