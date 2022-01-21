rule bowtie2_index:
    input:
        reference = config["reference"]
    params:
        basename = reference_dir_path / reference_basename
    output:
        index1=expand(reference_dir_path / "{basename}.{index}.bt2", basename = reference_basename, index = range(1, 5)),
        index2=expand(reference_dir_path / "{basename}.rev.{index}.bt2", basename = reference_basename, index = range(1, 3))
    log:
        std=log_dir_path / "bowtie2_index.log",
        cluster_log=cluster_log_dir_path / "bowtie2_index.cluster.log",
        cluster_err=cluster_log_dir_path / "bowtie2_index.cluster.err"
    benchmark:
        benchmark_dir_path / "/bowtie2_index.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bowtie2_threads"],
        time=config["bowtie2_time"],
        mem=config["bowtie2_mem_mb"]
    threads:
        config["bowtie2_threads"]
    shell:
        "bowtie2-build {input} {params.basename} > {log.std} 2>&1 || true"


rule samtools_bam_index:
    input:
        bam_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam",
        bam_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam"
    output:
        bai_raw=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam.bai",
        bai_clipped=clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam.bai"
    log:
        std=log_dir_path / "{sample_id}/bam_index.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/bam_index.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/bam_index.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/bam_index.benchmark.txt"
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


rule samtools_faidx:
    input:
        config["reference"]
    output:
        "%s.fai" % config["reference"]
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


rule picard_dict:
    input:
        config["reference"]
    output:
        f"{reference_dir_path}/{reference_basename}.dict"
    log:
        std=log_dir_path / "dict.log",
        cluster_log=cluster_log_dir_path / "dict.cluster.log",
        cluster_err=cluster_log_dir_path / "dict.cluster.err"
    benchmark:
        benchmark_dir_path / "dict.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["dict_threads"],
        time=config["dict_time"],
        mem=config["dict_mem_mb"],
    threads:
        config["dict_threads"]
    shell:
         "picard CreateSequenceDictionary R={input} O={output} > {log.std} 2>&1"