rule gatk_mutect2:
    input:
        reference=reference,
        dict=rules.picard_dict.output,
        fai=rules.samtools_faidx.output,
        sample=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam",
        index=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam.bai"
    output:
        vcf=varcall_gatk_dir_path / "{sample_id}/{sample_id}.mutect2.vcf.gz",
        stats=varcall_gatk_dir_path / "{sample_id}/{sample_id}.mutect2.vcf.gz.stats",
        tbi=varcall_gatk_dir_path / "{sample_id}/{sample_id}.mutect2.vcf.gz.tbi"
    params:
        options=config["gatk_mutect2_options"]
    log:
        std=log_dir_path / "{sample_id}/gatk_mutect2.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/gatk_mutect2.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/gatk_mutect2.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/gatk_mutect2.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["gatk_mutect2_threads"],
        mem=config["gatk_mutect2_mem_mb"],
        time=config["gatk_mutect2_time"]
    threads:
        config["gatk_mutect2_threads"]
    shell:
        "gatk --java-options '-Xmx{resources.mem}m' Mutect2 {params.options} -R {input.reference} -I {input.sample} "
        "-O {output.vcf} 1> {log.std} 2>&1 "


rule gatk_mergevcfs:
    input:
        expand(varcall_gatk_dir_path / "{sample_id}/{sample_id}.mutect2.vcf.gz", sample_id=config["sample_id"])
    output:
        varcall_gatk_dir_path / gatk_mergedvcf_filename
    log:
        std=log_dir_path / "gatk_mergevcfs.log",
        cluster_log=cluster_log_dir_path / "gatk_mergevcfs.cluster.log",
        cluster_err=cluster_log_dir_path / "gatk_mergevcfs.cluster.err"
    benchmark:
        benchmark_dir_path / "gatk_mergevcfs.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["gatk_mergevcfs_threads"],
        mem=config["gatk_mergevcfs_mem_mb"],
        time=config["gatk_mergevcfs_time"]
    threads:
        config["gatk_mergevcfs_threads"]
    shell:
        "for FILE in {input}; do echo $FILE >> tmp.list; done; "
        "gatk --java-options '-Xmx{resources.mem}m' MergeVcfs -I tmp.list -O {output} 1> {log.std} 2>&1; "
        "rm tmp.list "