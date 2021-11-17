rule pisces_somatic:
    input:
        ref_dir = reference_dir_path,
        sample = clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam",
        index = clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam.bai"
    output:
        directory(varcall_pisces_dir_path / "somatic/{sample_id}")
    params:
        pisces_tool_path=config["pisces_tool_path"],
        options=config["pisces_somatic_options"]
    log:
        std=log_dir_path / "{sample_id}/pisces_somatic.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/pisces_somatic.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/pisces_somatic.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/pisces_somatic.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["pisces_somatic_threads"],
        mem=config["pisces_somatic_mem_mb"],
        time=config["pisces_somatic_time"]
    threads:
        config["pisces_somatic_threads"]
    shell:
        "{params.pisces_tool_path}/Pisces -bam {input.sample} -g {input.ref_dir} {params.options} -OutFolder {output}"


rule pisces_germline:
    input:
        ref_dir = reference_dir_path,
        sample = clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam",
        index = clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam.bai"
    output:
        directory(varcall_pisces_dir_path / "germline/{sample_id}")
    params:
        pisces_tool_path = config["pisces_tool_path"],
        options = config["pisces_somatic_options"]
    log:
        std = log_dir_path / "{sample_id}/pisces_germline.log",
        cluster_log = cluster_log_dir_path / "{sample_id}/pisces_germline.cluster.log",
        cluster_err = cluster_log_dir_path / "{sample_id}/pisces_germline.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/pisces_germline.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus = config["pisces_germline_threads"],
        mem = config["pisces_germline_mem_mb"],
        time = config["pisces_germline_time"]
    threads:
        config["pisces_germline_threads"]
    shell:
        "{params.pisces_tool_path}/Pisces -bam {input.sample} -g {input.ref_dir} {params.options} -OutFolder {output}"