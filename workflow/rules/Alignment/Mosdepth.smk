rule mosdepth:
    input:
        bam=rules.bowtie2_map.output.bam,
        # bai=rules.index_bam.output
    output:
        alignment_dir_path / "{sample_id}/{sample_id}.coverage.per-base.bed.gz"
    params:
        min_mapping_quality=config["mosdepth_min_mapping_quality"],
        output_pefix=alignment_dir_path / "{sample_id}/{sample_id}.coverage"
    log:
        std=log_dir_path / "{sample_id}/mosdepth.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.mosdepth.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.mosdepth.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample_id}/mosdepth.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads:
        config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.output_pefix} {input.bam} > {log.std} 2>&1"