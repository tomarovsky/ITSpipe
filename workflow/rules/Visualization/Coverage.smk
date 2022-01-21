rule draw_coverage_raw:
    input:
        coverage_raw=raw_coverage_dir_path / "{sample_id}.sort.genomecov.tab.gz",
    output:
        coverage_raw_plot=raw_coverage_dir_path / "{sample_id}.plot.{ext}",
    params:
        tool=config["draw_coverage_tool"],
        output_raw_prefix=lambda wildcards, output: output["coverage_raw_plot"][:-9],
        extension=lambda wildcards, output: output["coverage_raw_plot"][-3:],
        xlabel=config["draw_coverage_xlabel"],
        ylabel=config["draw_coverage_ylabel"],
        title=config["draw_coverage_title"],
        width=config["draw_coverage_width"],
        height=config["draw_coverage_height"],
        markersize=config["draw_coverage_markersize"],
        ylogbase=config["draw_coverage_ylogbase"],
        type=config["draw_coverage_type"],
        grid=config["draw_coverage_grid"]
    log:
        std=log_dir_path / "{sample_id}/draw_coverage_raw.{ext}.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/draw_coverage_raw.{ext}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/draw_coverage_raw.{ext}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/draw_coverage_raw.{ext}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["visualization_threads"],
        time=config["visualization_time"],
        mem=config["visualization_mem_mb"],
    threads:
        config["visualization_threads"]
    shell:
        "python3 workflow/scripts/draw_coverage.py --input-file {input.coverage_raw} "
        "--output-prefix {params.output_raw_prefix} "
        "--tool '{params.tool}' "
        "--extensions '{params.extension}' "
        "--xlabel {params.xlabel} "
        "--ylabel {params.ylabel} "
        "--title '{params.title}' "
        "--width {params.width} "
        "--height {params.height} "
        "--markersize {params.markersize} "
        "--ylogbase {params.ylogbase} "
        "--type {params.type} "
        "--grid {params.grid} > {log.std} 2>&1; "


rule draw_coverage_clipped:
    input:
        coverage_clipped=clipped_coverage_dir_path / "{sample_id}.clipped.genomecov.tab.gz"
    output:
        coverage_clipped_plot=clipped_coverage_dir_path / "{sample_id}.clipped.trim.plot.{ext}"
    params:
        tool = config["draw_coverage_tool"],
        output_clipped_prefix=lambda wildcards, output: output["coverage_clipped_plot"][:-9],
        extension=lambda wildcards, output: output["coverage_clipped_plot"][-3:],
        xlabel=config["draw_coverage_xlabel"],
        ylabel=config["draw_coverage_ylabel"],
        title=config["draw_coverage_title"],
        width=config["draw_coverage_width"],
        height=config["draw_coverage_height"],
        markersize=config["draw_coverage_markersize"],
        ylogbase=config["draw_coverage_ylogbase"],
        type=config["draw_coverage_type"],
        grid=config["draw_coverage_grid"]
    log:
        std=log_dir_path / "{sample_id}/draw_coverage_clipped.{ext}.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/draw_coverage_clipped.{ext}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/draw_coverage_clipped.{ext}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/draw_coverage_clipped.{ext}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["visualization_threads"],
        time=config["visualization_time"],
        mem=config["visualization_mem_mb"],
    threads:
        config["visualization_threads"]
    shell:
        "python3 workflow/scripts/draw_coverage.py --input-file {input.coverage_clipped} "
        "--output-prefix {params.output_clipped_prefix} "
        "--tool '{params.tool}' "
        "--extensions '{params.extension}' "
        "--xlabel {params.xlabel} "
        "--ylabel {params.ylabel} "
        "--title '{params.title}' "
        "--width {params.width} "
        "--height {params.height} "
        "--markersize {params.markersize} "
        "--ylogbase {params.ylogbase} "
        "--type {params.type} "
        "--grid {params.grid} > {log.std} 2>&1; "