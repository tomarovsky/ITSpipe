rule ref_dict:
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
         "picard CreateSequenceDictionary -R {input} > {log.std} 2>&1"