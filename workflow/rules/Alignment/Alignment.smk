from Snakefile import raw_alignment_dir_path

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
        bowtie2_threads=config["bowtie2_threads"],
        tmp_prefix=lambda wildcards, output: output["sam"][:-4],
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
        config["bowtie2_threads"] # + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"]
    shell:
        "mkdir -p {output.outdir} ; "
        "bowtie2 -p {params.bowtie2_threads} -x {params.reference} -1 {input.forward_reads} -2 {input.reverse_reads} "
        "--rg-id '{wildcards.sample_id}' --rg 'ID:{wildcards.sample_id}' "
        "--rg 'SM:{wildcards.sample_id}' --rg 'PL:Illumina' --rg 'PU:x' --rg 'LB:x' -S {output.sam} 2> {log.bowtie2} "


rule sam_trimmer:
    input:
        rules.bowtie2_map.output.sam
    output:
        raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sam"
    params:
        pattern=config["sam_trimmer_pattern"],
        reference_start=config["sam_trimmer_reference_start"]
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
        "python3 workflow/scripts/SAM_trimmer.py -i {input} --pattern {params.pattern} --reference_start {params.reference_start} -o {output} > {log.std} 2>&1"


rule samtools_bam_improvements:
    input:
        sam=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sam"
    output:
        bam=raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam"
    params:
        fixmate_threads=config["fixmate_threads"],
        sort_threads=config["sort_threads"],
        markdup_threads=config["markdup_threads"],
        view_threads=config["view_threads"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"],
        view=config["samtools_view_options"],
        tmp_prefix=lambda wildcards, output: output["bam"][:-4],
    log:
        fixmate=log_dir_path / "{sample_id}/fixmate.log",
        sort=log_dir_path / "{sample_id}/sort.log",
        view=log_dir_path / "{sample_id}/view.log",
        markdup=log_dir_path / "{sample_id}/markdup.log",
        cluster_log=cluster_log_dir_path / "{sample_id}/samtools_bam_improvements.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}/samtools_bam_improvements.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/samtools_bam_improvements.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["sort_threads"] + config["fixmate_threads"] + config["view_threads"] + 1,  # + config["markdup_threads"],
        mem=config["per_thread_sort_mem"] * config["sort_threads"] * 1024 + config["fixmate_mem_mb"] + config["view_mem_mb"], # + config["markdup_mem_mb"],
        time=config["bowtie2_time"]
    threads:
        config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"]
    shell:
        # "samtools fixmate -@ {params.fixmate_threads} -m - - 2> {log.fixmate} | "
        "samtools sort -T {params.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} {input.sam} 2> {log.sort} | "
        # "samtools markdup -@ {params.markdup_threads} - 2> {log.markdup} | "
        "samtools view {params.view} -@ {params.view_threads} -o {output.bam} - 2> {log.view}; "