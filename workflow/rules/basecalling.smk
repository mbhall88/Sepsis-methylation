rule basecall:
    input:
        fast5_dir=get_fast5,
    output:
        seqsum=join(
            "results", module_name, rule_name, "{cond}_{rep}", "sequencing_summary.txt"
        ),
        fastq_dir=directory(join("results", module_name, rule_name, "{cond}_{rep}")),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=get_mem(config, rule_name),
    container:
        containers["guppy"]
    shell:
        "guppy_basecaller {params.opt} -i {input.fast5_dir} -s {output.fastq_dir} &> {log}"


rule_name = "merge_fastq"


rule merge_fastq:
    input:
        fastq_dir=rules.ont_guppy.output.fastq_dir,
    output:
        fastq=join("results", module_name, rule_name, "{cond}_{rep}.fastq"),
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
        "pyBioTools Fastq Filter {params.opt} -i {input.fastq_dir} -o {output.fastq} --verbose &> {log}"
