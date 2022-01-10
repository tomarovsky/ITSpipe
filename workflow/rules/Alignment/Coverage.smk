rule mosdepth:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.bam",
        bai_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.bam.bai",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam",
        bai_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam.bai"
    output:
        coverage_raw=raw_coverage_dir_path / "{sample_id}.coverage.per-base.bed.gz",
        coverage_clipped=clipped_coverage_dir_path / "{sample_id}.clipped.coverage.per-base.bed.gz"
    params:
        min_mapping_quality=config["mosdepth_min_mapping_quality"],
        output_raw_pefix=lambda wildcards, output: output["coverage_raw"][:-16],
        output_clipped_pefix=lambda wildcards, output: output["coverage_clipped"][:-16],
    log:
        std=log_dir_path / "{sample_id}/mosdepth.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/mosdepth.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/mosdepth.cluster.err"
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
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.output_raw_pefix} {input.bam_raw} > {log.std} 2>&1; "
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.output_clipped_pefix} {input.bam_clipped} > {log.std} 2>&1; "


rule genomecov:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam"
    output:
        coverage_raw=raw_coverage_dir_path / "{sample_id}.trim.sort.genomecov.tab.gz",
        coverage_clipped=clipped_coverage_dir_path / "{sample_id}.clipped.trim.genomecov.tab.gz"
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