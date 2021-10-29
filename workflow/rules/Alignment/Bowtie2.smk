# localrules: bowtie2_index
# ruleorder: bowtie2_index > bowtie2_map


rule bowtie2_map:
    input:
        forward_reads=rules.trimmomatic.output.pe_forward,
        reverse_reads=rules.trimmomatic.output.pe_reverse,
        reference=config["reference"]
    output:
        bam=alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam"
    params:
        fixmate_threads=config["fixmate_threads"],
        sort_threads=config["sort_threads"],
        markdup_threads=config["markdup_threads"],
        bowtie2_threads=config["bowtie2_threads"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"],
        tmp_prefix=alignment_dir_path / "{sample_id}/{sample_id}",
        outdir = alignment_dir_path / "{sample_id}"
    log:
        bowtie2=log_dir_path / "{sample_id}/bowtie2.log",
        fixmate=log_dir_path / "{sample_id}/fixmate.log",
        sort=log_dir_path / "{sample_id}/sort.log",
        markdup=log_dir_path / "{sample_id}/markdup.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.bowtie2_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.bowtie2_map.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/bowtie2_map.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bowtie2_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"] + 1,
        mem=config["per_thread_sort_mem"] * config["sort_threads"] * 1024 + config["bowtie2_mem_mb"] + config["fixmate_mem_mb"] + config["markdup_mem_mb"],
        time=config["bowtie2_time"]
    threads:
        config["bowtie2_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"]
    shell:
        """
        mkdir -r {params.outdir}; \
        bowtie2 -p {params.bowtie2_threads} {input.reference} -1 <(gunzip -c {input.forward_reads}) -2 <(gunzip -c {input.reverse_reads}) \
        --rg \'@RG\\tID:{wildcards.sample_id}\\tPU:x\\tSM:{wildcards.sample_id}\\tPL:Illumina\\tLB:x\' 2> {log.bowtie2} | \
        samtools fixmate -@ {params.fixmate_threads} -m - - 2> {log.fixmate} | \
        samtools sort -T {params.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} 2> {log.sort} | \
        samtools markdup -@ {params.markdup_threads} - {output.bam} 2> {log.markdup}
        """
#
# rule bowtie2_index:
#     input:
#         reference = config["reference"]
#     params:
#         basename = reference_dir_path / reference_basename
#     output:
#         expand(reference_dir_path / "{basename}.{index}.bt2", basename = reference_basename, index = range(1, 5)),
#         expand(reference_dir_path / "{basename}.rev.{index}.bt2", basename = reference_basename, index = range(1, 3))
#     log:
#         std=log_dir_path / "bowtie2_index.log",
#         cluster_log=cluster_log_dir_path / "bowtie2_index.cluster.log",
#         cluster_err=cluster_log_dir_path / "bowtie2_index.cluster.err"
#     benchmark:
#         benchmark_dir_path / "/bowtie2_index.benchmark.txt"
#     conda:
#         "../../../%s" % config["conda_config"]
#     resources:
#         cpus=config["bowtie2_threads"],
#         time=config["bowtie2_time"],
#         mem=config["bowtie2_mem_mb"]
#     threads:
#         config["bowtie2_threads"]
#     shell:
#         "bowtie2-build {input} {params.basename}"