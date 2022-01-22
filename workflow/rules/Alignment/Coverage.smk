rule genomecov:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.sort.bam",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam"
    output:
        coverage_raw=raw_coverage_dir_path / "{sample_id}.sort.genomecov.tab.gz",
        coverage_clipped=clipped_coverage_dir_path / "{sample_id}.clipped.genomecov.tab.gz"
    params:
        options=config["bedtools_genomecov_options"],
    log:
        std=log_dir_path / "{sample_id}/genomecov.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/genomecov.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/genomecov.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/genomecov.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bedtools_genomecov_threads"],
        time=config["bedtools_genomecov_time"],
        mem=config["bedtools_genomecov_mem_mb"],
    threads:
        config["bedtools_genomecov_threads"]
    shell:
        "bedtools genomecov -ibam {input.bam_raw} {params.options} | gzip > {output.coverage_raw} 2> {log.std}; "
        "bedtools genomecov -ibam {input.bam_clipped} {params.options} | gzip > {output.coverage_clipped} 2> {log.std}; "