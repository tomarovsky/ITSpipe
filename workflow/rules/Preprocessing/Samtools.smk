rule samtools_bam_index:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.bam"
    output:
        bai_raw=temp(raw_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam.bai"),
        bai_clipped=temp(clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.bam.bai")
    log:
        std=log_dir_path / "{sample_id}.bam_index.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.bam_index.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.bam_index.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.bam_index.benchmark.txt"
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


rule samtools_view:
    input:
        bam=clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.bam"
    output:
        view_bam=clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam"
    params:
        options=config["samtools_view_options"]
    log:
        std=log_dir_path / "{sample_id}.view.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.view.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.view.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.view.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["view_threads"],
        time=config["view_time"],
        mem=config["view_mem_mb"],
    threads:
        config["view_threads"]
    shell:
        "samtools view {params.options} -@ {threads} -o {output.view_bam} {input.bam} > {log.std} 2>&1; "


rule samtools_view_bam_index:
    input:
        view_bam=clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam"
    output:
        view_bam_bai=temp(raw_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam.bai")
    log:
        std=log_dir_path / "{sample_id}.view_bam_index.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.view_bam_index.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.view_bam_index.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.view_bam_index.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["index_threads"],
        time=config["index_time"],
        mem=config["index_mem_mb"],
    threads:
        config["index_threads"]
    shell:
        "samtools index -@ {threads} {input.view_bam} > {log.std} 2>&1; "


rule samtools_faidx:
    input:
        reference
    output:
        temp(f"{reference}.fai")
    log:
        std=log_dir_path / "faidx.log",
        cluster_log=cluster_log_dir_path / "faidx.cluster.log",
        cluster_err=cluster_log_dir_path / "faidx.cluster.err"
    benchmark:
        benchmark_dir_path / "faidx.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["index_threads"],
        time=config["index_time"],
        mem=config["index_mem_mb"],
    threads:
        config["index_threads"]
    shell:
        "samtools faidx {input} > {log.std} 2>&1; "