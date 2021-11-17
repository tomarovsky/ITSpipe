rule bamutil_clipoverlap:
    input:
        rules.bowtie2_map.output.bam
    output:
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam"
    params:
        poolsize=config["poolsize"]
    log:
        std=log_dir_path / "{sample_id}/bamutil_clipoverlap.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/bamutil_clipoverlap.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/bamutil_clipoverlap.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample_id}/bamutil_clipoverlap.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bamutil_clipoverlap_threads"],
        time=config["bamutil_clipoverlap_time"],
        mem=config["bamutil_clipoverlap_mem_mb"],
    threads:
        config["bamutil_clipoverlap_threads"]
    shell:
        "bam clipOverlap --in {input} --out {output} --poolSize {params.poolsize} > {log.std} 2>&1"