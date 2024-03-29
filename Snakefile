from pathlib import Path


#---- setup config ----
configfile: "config/default.yaml"

#---- setup paths ----
cluster_log_dir_path = Path(config["cluster_log_dir"])
log_dir_path = Path(config["log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])
samples_dir_path = Path(config["samples_dir"])
output_dir_path = Path(config["output_dir"])

filtered_reads_dir_path = output_dir_path / config["filtered_read_dir"]
raw_alignment_dir_path = output_dir_path / config["raw_alignment_dir"]
raw_coverage_dir_path = output_dir_path / config["raw_coverage_dir"]
clipped_alignment_dir_path = output_dir_path / config["clipped_alignment_dir"]
clipped_coverage_dir_path = output_dir_path / config["clipped_coverage_dir"]
varcall_bcftools_dir_path = output_dir_path / config["varcall_bcftools_dir"]
varcall_gatk_dir_path = output_dir_path / config["varcall_gatk_dir"]
varcall_pisces_dir_path = output_dir_path / config["varcall_pisces_dir"]

#---- setup filenames ----
reference = Path(config["reference"])
reference_dir_path = reference.parents[0]
reference_filename = reference.name
reference_basename = reference.stem
gatk_merged_vcf_filename = "%s.%s.vcf.gz" % (reference_basename, config["gatk_merged_vcf_prefix"])
pisces_somatic_merged_vcf_filename = "%s.%s.vcf.gz" % (reference_basename, config["pisces_somatic_merged_vcf_prefix"])
pisces_germline_merged_vcf_filename = "%s.%s.vcf.gz" % (reference_basename, config["pisces_germline_merged_vcf_prefix"])


if "sample_id" not in config:
    config["sample_id"] = [d.name for d in samples_dir_path.iterdir() if d.is_dir()]
if "draw_coverage_plot_extensions" in config:
    config["draw_coverage_plot_extensions"] = [ext for ext in config["draw_coverage_plot_extensions"].strip().split(",")]


localrules: all

rule all:
    input:
        # Trimmomatic:
        expand(filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_1.fastq.gz", sample_id=config["sample_id"]),
        expand(filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_1.se.fastq.gz", sample_id=config["sample_id"]),
        expand(filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_2.fastq.gz", sample_id=config["sample_id"]),
        expand(filtered_reads_dir_path / "{sample_id}/{sample_id}.trimmed_2.se.fastq.gz", sample_id=config["sample_id"]),

        # Bowtie2:
        expand(raw_alignment_dir_path / "{sample_id}/{sample_id}.sam", sample_id=config["sample_id"]),

        # BAM_trimmer:
        expand(raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sam", sample_id=config["sample_id"]),

        # Samtools BAM improvements:
        expand(raw_alignment_dir_path / "{sample_id}/{sample_id}.sort.bam", sample_id=config["sample_id"]),
        expand(raw_alignment_dir_path / "{sample_id}/{sample_id}.trim.sort.bam", sample_id=config["sample_id"]),

        # Bamutil:
        expand(clipped_alignment_dir_path / "{sample_id}/{sample_id}.clipped.bam", sample_id=config["sample_id"]),

        # Genomecov:
        expand(raw_coverage_dir_path / "{sample_id}.sort.genomecov.tab.gz", sample_id=config["sample_id"]),
        expand(clipped_coverage_dir_path / "{sample_id}.clipped.genomecov.tab.gz", sample_id=config["sample_id"]),

        # Coverage visualization:
        expand(raw_coverage_dir_path / "{sample_id}.plot.{ext}", sample_id=config["sample_id"], ext=config["draw_coverage_plot_extensions"]),
        expand(clipped_coverage_dir_path / "{sample_id}.clipped.trim.plot.{ext}", sample_id=config["sample_id"], ext=config["draw_coverage_plot_extensions"]),

        # Variant calling:
        # Bcftools:
        expand(varcall_bcftools_dir_path / "{reference_basename}.vcf.gz", reference_basename = reference_basename),
        expand(varcall_bcftools_dir_path / "{reference_basename}.mpileup.vcf.gz", reference_basename = reference_basename),
        expand(varcall_bcftools_dir_path / "{reference_basename}.mpileup.filt.vcf.gz", reference_basename = reference_basename),
        # GATK:
        expand(varcall_gatk_dir_path / "{sample_id}/{sample_id}.mutect2.vcf.gz", sample_id=config["sample_id"]),
        varcall_gatk_dir_path / gatk_merged_vcf_filename,
        # Pisces:
        expand(varcall_pisces_dir_path / "somatic/{sample_id}/{sample_id}.clipped.vcf.gz", sample_id=config["sample_id"]),
        expand(varcall_pisces_dir_path / "germline/{sample_id}/{sample_id}.clipped.vcf.gz", sample_id=config["sample_id"]),
        varcall_pisces_dir_path / pisces_somatic_merged_vcf_filename,
        varcall_pisces_dir_path / pisces_germline_merged_vcf_filename


#---- load rules ----
include: "workflow/rules/QCFiltering/Trimmomatic.smk"
include: "workflow/rules/Preprocessing/Indexes.smk"
include: "workflow/rules/Alignment/Alignment.smk"
include: "workflow/rules/Alignment/Coverage.smk"
include: "workflow/rules/QCFiltering/Bamutil.smk"
include: "workflow/rules/Visualization/Coverage.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"
include: "workflow/rules/VariantCall/Gatk.smk"
include: "workflow/rules/VariantCall/Pisces.smk"