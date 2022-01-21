localrules: bowtie2_index
ruleorder: bowtie2_index > bowtie2_map


rule bowtie2_map:
    input:
        forward_reads=rules.trimmomatic.output.pe_forward,
        reverse_reads=rules.trimmomatic.output.pe_reverse,
        reference=config["reference"],
        index1=expand(reference_dir_path / "{basename}.{index}.bt2", basename = reference_basename, index = range(1, 5)),
        index2=expand(reference_dir_path / "{basename}.rev.{index}.bt2", basename = reference_basename, index = range(1, 3))
    output:
        outdir=directory(raw_alignment_dir_path/"{sample_id}"),
        sam=raw_alignment_dir_path / "{sample_id}/{sample_id}.sam"
    params:
        reference=reference_dir_path / reference_basename
    log:
        bowtie2=log_dir_path / "{sample_id}/bowtie2.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/bowtie2_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/bowtie2_map.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/bowtie2_map.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bowtie2_threads"],
        mem=config["bowtie2_mem_mb"],
        time=config["bowtie2_time"]
    threads:
        config["bowtie2_threads"]
    shell:
        "mkdir -p {output.outdir} ; "
        "bowtie2 -p {threads} -x {params.reference} -1 {input.forward_reads} -2 {input.reverse_reads} "
        "--rg-id '{wildcards.sample_id}' --rg 'ID:{wildcards.sample_id}' "
        "--rg 'SM:{wildcards.sample_id}' --rg 'PL:Illumina' --rg 'PU:x' --rg 'LB:x' -S {output.sam} 2> {log.bowtie2} "


rule SAM_trimmer:
    input:
        rules.bowtie2_map.output.sam
    output:
        raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sam"
    params:
        pattern=config["sam_trimmer_pattern"],
        pos=config["sam_trimmer_pos"]
    log:
        std=log_dir_path / "{sample_id}/sam_trimmer.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/sam_trimmer.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/sam_trimmer.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample_id}/sam_trimmer.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["sam_trimmer_threads"],
        time=config["sam_trimmer_time"],
        mem=config["sam_trimmer_mem_mb"],
    threads:
        config["sam_trimmer_threads"]
    shell:
        "python3 workflow/scripts/SAM_trimmer.py -i {input} --pattern {params.pattern} --pos {params.pos} -o {output} > {log.std} 2>&1"


rule samtools_bam_improvements:
    input:
        sam=raw_alignment_dir_path / "{sample_id}/{sample_id}.sam",
        trim_sam=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sam"
    output:
        bam=raw_alignment_dir_path / "{sample_id}/{sample_id}.sort.bam",
        trim_bam=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam",
    params:
        sort_threads=config["sort_threads"],
        view_threads=config["view_threads"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"],
        view=config["samtools_view_options"],
        tmp_prefix=lambda wildcards, output: output["bam"][:-4],
    log:
        sort=log_dir_path / "{sample_id}/samtools_bam_improvements.sort.log",
        view=log_dir_path / "{sample_id}/samtools_bam_improvements.view.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/samtools_bam_improvements.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/samtools_bam_improvements.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/samtools_bam_improvements.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["sort_threads"] + config["view_threads"],
        mem=config["per_thread_sort_mem"] * config["sort_threads"] * 1024 + config["view_mem_mb"],
        time=config["samtools_bam_improvements_time"]
    threads:
        config["sort_threads"] + config["view_threads"]
    shell:
        "samtools sort -T {params.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} {input.sam} 2> {log.sort} | "
        "samtools view {params.view} -@ {params.view_threads} -o {output.bam} - 2> {log.view}; "
        "samtools sort -T {params.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} {input.trim_sam} 2> {log.sort} | "
        "samtools view {params.view} -@ {params.view_threads} -o {output.trim_bam} - 2> {log.view}; "