rule bam_trimmer:
    input:
        rules.bamutil_clipoverlap.output.bam_clipped
    output:
        clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.trim.bam"
    params:
        pattern=config["bam_trimmer_pattern"],
        reference_start=config["bam_trimmer_reference_start"]
    log:
        std=log_dir_path / "{sample_id}/bam_trimmer.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/bam_trimmer.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/bam_trimmer.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample_id}/bam_trimmer.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bam_trimmer_threads"],
        time=config["bam_trimmer_time"],
        mem=config["bam_trimmer_mem_mb"],
    threads:
        config["bam_trimmer_threads"]
    shell:
        "python3 workflow/scripts/BAM_trimmer.py -i {input} --pattern {params.pattern} --reference_start {params.reference_start} -o {output} > {log.std} 2>&1"