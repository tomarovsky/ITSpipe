rule pisces_somatic:
    input:
        ref_dir = reference_dir_path,
        sample = clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam",
        index = clipped_alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.clipped.view.bam.bai"
    output:
        directory(varcall_pisces_dir_path / "{sample_id}")
    params:
        options=config["pisces_somatic_options"]
    log:
        std=log_dir_path / "{sample_id}.pisces_somatic.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.pisces_somatic.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.pisces_somatic.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}.pisces_somatic.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["pisces_somatic_threads"],
        mem=config["pisces_somatic_mem_mb"],
        time=config["pisces_somatic_time"]
    threads:
        config["pisces_somatic_threads"]
    shell:
        "dotnet Pisces.dll -bam {input.sample} -g {input.ref_dir} {params.options} -OutFolder {output}"





# dotnet Pisces.dll -bam /my/path/to/TestData/example_S1.bam -g /my/path/to/WholeGenomeFasta
# Somatic: -bam {Bam} -CallMNVs false -g {genome folder} -gVCF false -i {intervalfile} -OutFolder {outfolder}
# Germline: -bam {Bam} -CallMNVs false -crushvcf true -g {genomefolder} -gVCF false -i {intervalfile} -ploidy diploid -OutFolder {outfolder}