localrules: trimmomatic


rule trimmomatic:
    input:
        samples_dir_path / "{sample_id}/{sample_id}_1.fastq.gz",
        samples_dir_path / "{sample_id}/{sample_id}_2.fastq.gz"
    output:
        pe_forward=filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_1.fastq.gz",
        se_forward=filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_1.se.fastq.gz",
        pe_reverse=filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_2.fastq.gz",
        se_reverse=filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_2.se.fastq.gz"
    params:
        adapters=config["adapters"],
        illumina_clip="2:30:10:1",
        window_size=8,
        window_quality=20,
        minlength=50
    log:
        std=log_dir_path / "{sample_id}/trimmomatic.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.trimmomatic.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.trimmomatic.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/trimmomatic.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["trimmomatic_threads"],
        time=config["trimmomatic_time"],
        mem=config["trimmomatic_mem_mb"]
    threads:
        config["trimmomatic_threads"]
    shell:
         "trimmomatic PE -threads {threads} -phred33 {input} {output} "
         "ILLUMINACLIP:{params.adapters}:{params.illumina_clip} "
         "SLIDINGWINDOW:{params.window_size}:{params.window_quality} "
         "MINLEN:{params.minlength} > {log.std} 2>&1"