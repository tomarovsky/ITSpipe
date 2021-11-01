rule index_bam:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.bam"
    output:
        bai_raw=temp(alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam.bai"),
        bai_clipped=temp(alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.bam.bai")
    log:
        std=log_dir_path / "{sample_id}/index.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.index.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.index.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/index.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["index_threads"],
        time=config["index_time"],
        mem=config["index_mem_mb"],
    threads:
        config["index_threads"]
    shell:
        "samtools index -@ {threads} {input.bam_raw} > {log.std} 2>&1; "
        "samtools index -@ {threads} {input.bam_clipped} > {log.std} 2>&1; "