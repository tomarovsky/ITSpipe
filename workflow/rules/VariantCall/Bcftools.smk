rule bcftools_mpileup:
    input:
        ref=reference,
        samples=expand(clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam", sample_id=config["sample_id"]),
        indexes=expand(clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam.bai", sample_id=config["sample_id"])
    output:
        varcall_dir_path / "{reference_basename}.mpileup.vcf.gz"
    log:
        mpileup=log_dir_path / "{reference_basename}.mpileup.log",
        call=log_dir_path / "{reference_basename}.call.log",
        cluster_log=cluster_log_dir_path / "{reference_basename}.bowtie2_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{reference_basename}.bowtie2_map.cluster.err"
    benchmark:
        benchmark_dir_path / "{reference_basename}.bcftools_mpileup.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_mpileup_threads"],
        mem=config["bcftools_mpileup_mem_mb"],
        time=config["bcftools_mpileup_time"]
    threads:
        config["bcftools_mpileup_threads"] #-d 250 -q 30 -Q 30
    shell:
        "bcftools mpileup --threads {threads} --adjust-MQ 50 -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR -Ou -f {input.ref} {input.samples} | "
        "bcftools call -m -O z -v -f GQ,GP > {output}; "