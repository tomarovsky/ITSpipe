rule index_bam:
    input:
        rules.bowtie2_map.output.bam
    output:
        temp(alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam.bai")
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
        "samtools index -@ {threads} {input} > {log.std} 2>&1"