rule bcftools_mpileup:
    input:
        reference=reference,
        samples=expand(clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam", sample_id=config["sample_id"]),
        indexes=expand(clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam.bai", sample_id=config["sample_id"])
    output:
        pipe(varcall_bcftools_mpileup_dir_path / "{reference_basename}.mpileup.bcf")
    params:
        adjustMQ=50,
        annotate_mpileup=config["bcftools_mpileup_annotate"],
        # annotate_call=config["bcftools_call_annotate"],
        max_depth=config["bcftools_mpileup_max_depth"],
        min_MQ=config["bcftools_mpileup_min_MQ"],
        min_BQ=config["bcftools_mpileup_min_BQ"]
    # log:
    #     mpileup=log_dir_path / "{reference_basename}.bcftools_mpileup.log",
    #     call=log_dir_path / "{reference_basename}.bcftools_call.log",
    #     cluster_log=cluster_log_dir_path / "{reference_basename}.bcftools_mpileup.cluster.log",
    #     cluster_err=cluster_log_dir_path / "{reference_basename}.bcftools_mpileup.cluster.err"
    # benchmark:
    #     benchmark_dir_path / "{reference_basename}.bcftools_mpileup.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_mpileup_threads"],
        mem=config["bcftools_mpileup_mem_mb"],
        time=config["bcftools_mpileup_time"]
    threads:
        config["bcftools_mpileup_threads"]
    shell:
        "bcftools mpileup --threads {threads} -d {params.max_depth} -q {params.min_MQ} -Q {params.min_BQ} "
        "--adjust-MQ {params.adjustMQ} --annotate {params.annotate_mpileup} -Ou -f {input.reference} {input.samples} -o {output}"


rule bcftools_call:
    input:
        rules.bcftools_mpileup.output
    output:
        varcall_bcftools_mpileup_dir_path / "{reference_basename}.mpileup.vcf.gz"
    params:
        annotate_call=config["bcftools_call_annotate"],
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_mpileup_threads"],
        mem=config["bcftools_mpileup_mem_mb"],
        time=config["bcftools_mpileup_time"]
    threads:
        config["bcftools_mpileup_threads"]
    shell:
        "bcftools call -Oz -mv --annotate {params.annotate_call} -o {output} {input}"


rule bcftools_filter:
    input:
        rules.bcftools_call.output
    output:
        varcall_bcftools_mpileup_dir_path / "{reference_basename}.mpileup.filt.vcf.gz"
    params:
        soft_filter=config["bcftools_filter_soft_filter"],
        exclude=config["bcftools_filter_exclude"],
    log:
        std=log_dir_path / "{reference_basename}.bcftools_filter.log",
        cluster_log=cluster_log_dir_path / "{reference_basename}.bcftools_filter.cluster.log",
        cluster_err=cluster_log_dir_path / "{reference_basename}.bcftools_filter.cluster.err"
    benchmark:
        benchmark_dir_path / "{reference_basename}.bcftools_filter.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_threads"],
        mem=config["bcftools_filter_mem_mb"],
        time=config["bcftools_filter_time"]
    threads:
        config["bcftools_filter_threads"]
    shell:
        "bcftools filter -Oz -s {params.soft_filter} --exclude '{params.exclude}' {input} > {output} 2> {log.std}; "
