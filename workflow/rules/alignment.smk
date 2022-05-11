aln_dir = results / "alignment"


rule minimap2_index:
    input:
        fasta=rules.index_transcriptome.output.fasta,
    output:
        idx=data_dir / "reference/transcriptome_reference.mmi",
    log:
        logs_dir / "minimap2_index.log",
    params:
        opt="-x map-ont",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * GB,
        time="1h",
    conda:
        str(envs_dir / "aln_tools.yaml")
    shell:
        "minimap2 -t {threads} {params.opt} -d {output.idx} {input.fasta} 2> {log}"


rule minimap2_align:
    input:
        idx=rules.minimap2_index.output.idx,
        fastq=rules.merge_fastq.output.fastq,
    output:
        bam=temp(aln_dir / "{sample}.bam"),
        bam_index=temp(aln_dir / "{sample}.bam.bai"),
    log:
        logs_dir / "minimap2_align/{sample}.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
        time="2h",
    params:
        opt="-aL -x map-ont",
    conda:
        str(envs_dir / "aln_tools.yaml")
    shell:
        """
        (minimap2 -t {threads} {params.opt} {input.idx} {input.fastq} | \
            samtools sort -@ {threads} -o {output.bam}) 2> {log}
        samtools index {output.bam} 2>> {log}
        """


prefilter_dir = aln_dir / "prefilter"


rule alignmemt_prefilter:
    input:
        bam=rules.minimap2_align.output.bam,
        bam_index=rules.minimap2_align.output.bam_index,
    output:
        bam=temp(prefilter_dir / "{sample}.bam"),
        bam_index=temp(prefilter_dir / "{sample}.bam.bai"),
        reads_index=temp(prefilter_dir / "{sample}.bam.idx.gz"),
    resources:
        time="1h",
    log:
        logs_dir / "alignment_prefilter/{sample}.log",
    params:
        opt=" ".join(
            [
                "--skip_unmapped",
                "--skip_secondary",
                "--skip_supplementary",
                "--index_reads",
                "--orientation '+'",
                "--min_read_len 100",
                "--min_align_len 100",
                "--min_mapq 10",
                "--min_freq_identity 0.8",
            ]
        ),
    container:
        containers["pybiotools"]
    shell:
        "pyBioTools Alignment Filter {params.opt} -i {input.bam} -o {output.bam} --verbose &> {log}"


rule min_ref_coverage:
    input:
        reads_index_list=expand(
            prefilter_dir / "{sample}.bam.idx.gz",
            sample=TESTS,
        ),
    output:
        ref_list=aln_dir / "valid_references_list.txt",
    log:
        logs_dir / "min_ref_coverage.log",
    resources:
        time="1h",
    params:
        opt=dict(min_cov=30),
    container:
        containers["metacompore_python"]
    script:
        str(scripts_dir / "min_ref_coverage.py")


postfilter_dir = aln_dir / "postfilter"


rule alignmemt_postfilter:
    input:
        bam=rules.alignmemt_prefilter.output.bam,
        bam_index=rules.alignmemt_prefilter.output.bam_index,
        ref_list=rules.min_ref_coverage.output.ref_list,
    output:
        bam=postfilter_dir / "{sample}.bam",
        bam_index=postfilter_dir / "{sample}.bam.bai",
        reads_index=postfilter_dir / "{sample}.bam.idx.gz",
        selected_reads_fn=postfilter_dir / "{sample}_selected_read_ids.txt",
    log:
        logs_dir / "alignment_postfilter/{sample}.log",
    params:
        opt="--index_reads",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * GB,
        time="1h",
    container:
        containers["pybiotools"]
    shell:
        """
        pyBioTools Alignment Filter {params.opt} -i {input.bam} \
            --select_ref_fn {input.ref_list} -o {output.bam} \
            -l {output.selected_reads_fn} --verbose &> {log}
        """
