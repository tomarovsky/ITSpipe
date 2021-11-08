rule draw_coverage_raw:
    input:
        coverage_raw=raw_coverage_dir_path / "{sample_id}.coverage.per-base.bed.gz"
    output:
        coverage_raw_png=raw_coverage_dir_path / "{sample_id}.plot.png",
        coverage_raw_svg=raw_coverage_dir_path / "{sample_id}.plot.svg"
    params:
        output_raw_prefix=lambda wildcards, output: output["coverage_raw_png"][:-9],
        start_column_index=1,
        stop_column_index=2,
        coverage_column_index=3,
        min_x=None,
        max_x=None,
        min_y=None,
        max_y=None,
        extensions=config["plot_extensions"],
        xlabel=None,
        ylabel=None,
        title="{sample_id} coverage plot",
        width=6,
        height=6,
        markersize=2,
        type="plot",
        grid=False,
        close_plot=True
    log:
        std=log_dir_path / "{sample_id}.draw_coverage_raw.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.draw_coverage_raw.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.draw_coverage_raw.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.draw_coverage_raw.benchmark.txt"
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
        "--start_column_index {params.start_column_index} "
        "--stop_column_index {params.stop_column_index} "
        "--coverage_column_index {params.coverage_column_index} "
        "--extensions {params.extensions} "
        "--min_x {params.min_x} "
        "--max_x {params.max_x} "
        "--min_y {params.min_y} "
        "--max_y {params.max_y} "
        "--xlabel {params.xlabel} "
        "--ylabel {params.ylabel} "
        "--title {params.title} "
        "--width {params.width} "
        "--height {params.height} "
        "--markersize {params.markersize} "
        "--type {params.type} "
        "--grid {params.grid} "
        "--close_plot {params.close_plot} > {log.std} 2>&1; "


rule draw_coverage_clipped:
    input:
        coverage_clipped=clipped_coverage_dir_path / "{sample_id}.clipped.coverage.per-base.bed.gz"
    output:
        coverage_clipped_png=clipped_coverage_dir_path / "{sample_id}.clipped.plot.png",
        coverage_clipped_svg=clipped_coverage_dir_path / "{sample_id}.clipped.plot.svg"
    params:
        output_clipped_prefix=lambda wildcards, output: output["coverage_clipped_png"][:-9],
        start_column_index=1,
        stop_column_index=2,
        coverage_column_index=3,
        min_x=None,
        max_x=None,
        min_y=None,
        max_y=None,
        extensions=config["plot_extensions"],
        xlabel=None,
        ylabel=None,
        title="{sample_id} coverage plot",
        width=6,
        height=6,
        markersize=2,
        type="plot",
        grid=False,
        close_plot=True
    log:
        std=log_dir_path / "{sample_id}.draw_coverage_clipped.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.draw_coverage_clipped.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.draw_coverage_clipped.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.draw_coverage_clipped.benchmark.txt"
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
        "--start_column_index {params.start_column_index} "
        "--stop_column_index {params.stop_column_index} "
        "--coverage_column_index {params.coverage_column_index} "
        "--extensions {params.extensions} "
        "--min_x {params.min_x} "
        "--max_x {params.max_x} "
        "--min_y {params.min_y} "
        "--max_y {params.max_y} "
        "--xlabel {params.xlabel} "
        "--ylabel {params.ylabel} "
        "--title {params.title} "
        "--width {params.width} "
        "--height {params.height} "
        "--markersize {params.markersize} "
        "--type {params.type} "
        "--grid {params.grid} "
        "--close_plot {params.close_plot} > {log.std} 2>&1; "