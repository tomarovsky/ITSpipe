# ITSpipe  
  
Pipeline for the analysis of ITS sequences from the ribosomal cluster.

The pipeline includes:
- FASTQ filtering using Trimmomatic  
- Read alignments using Bowtie2
- Coverage using Bedtools
- Trimming and filtering alignments usig Bamutil and BAM_trimmer (custom Python script for cropping areas with high coverage in BAM file)
- Сoverage visualization (before and after trimming) using Matplotlib
- Variant calling using Gatk, Bcftools and Pisces 