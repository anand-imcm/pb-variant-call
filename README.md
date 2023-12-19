# Variant Calling and Annotation using Hi-Fi Reads

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/anand-imcm/pb-variant-call/publish.yml)&nbsp;&nbsp;
![GitHub release (with filter)](https://img.shields.io/github/v/release/anand-imcm/pb-variant-call)&nbsp;&nbsp;
[![Open](https://img.shields.io/badge/Open-Dockstore-blue)](https://dockstore.org/workflows/github.com/anand-imcm/pb-variant-call:main?tab=info)

This repository contains a WDL-based workflow for variant calling and annotation using Hi-Fi reads. The workflow includes several steps such as alignment, variant calling, VCF filtering, VCF normalization, variant phasing, variant annotation, and structural variant calling.

> To import the workflow into your Terra workspace, click on the Dockstore badge, and select 'Terra' from the 'Launch with' widget on the Dockstore workflow page.

## Workflow Steps

- **Alignment**: The HiFi reads are aligned to a reference genome using `pbmm2`. The output is a BAM file that contains the alignments.

- **Variant Calling**: Variants are called from the alignments using `DeepVariant`. The output is a VCF file that contains the called variants.

- **VCF Filtering**: The variants in the VCF file are filtered using `bcftools -f PASS` to include only variants that have passed all filters.

- **VCF Normalization**: The called variants are normalized using `bcftools norm`. This step ensures that all variants are represented in a standard way.

- **Variant Phasing**: The "PASS" variants are phased using `whatshap phase`. The phasing is encoded in the "FORMAT" column of the VCF.

- **Variant Annotation**: The "PASS" variants are annotated using `VEP` tool. The output is a VCF file that contains all the additional annotation in the "INFO" column.

- **Structural Variant Calling**: The structural variants are called from the alignments using `pbsv`. The output is a VCF
  file that contains the called structural variants.

- **Summary**: Custom scripts are used to generate coverage depth plot and a summary report.

## Inputs

The main inputs to the workflow are:

- **required**

  - `reads_fastq_gz` : Input PacBio HiFi reads in .fastq.gz format.
  - `prefix` : Sample name. This will be used as prefix for all the output files.
  - `genome_ref` : Human reference genome .fasta file. The version being used is GRCh38 release110 ([source](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/)).
  - `genome_index_pbmm` : Reference index generated through pbmm2 in .mmi format.
  - `vep_cache` : VEP cache in .zip format. This cache is required for the `VEP` tool to function. The version being used is Ensembl GRCh38 release v110 ([source](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)).
  - `target_bed` : "Coordinates for the target (amplified) regions. (0-based bed file). Ensembl/Gencode gene model is currently being used."

- **optional**

  - `vep_version` : The version of the VEP tool to use. Default value: `release_110.1`. This should be compatible with the `vep_cache` version.
  - `deepvariant_num_shards` : The number of shards to use when running DeepVariant. Default value : `12`.
  - `deepvariant_version` : The version of the DeepVariant tool to use. Default value: `1.5.0`.

## Outputs

The main outputs from the workflow are categorized based on the steps of the workflow:

- **Alignment**

  - `raw_hifi_to_reference_alignment_log`: Log file for the alignment step.
  - `raw_hifi_reads_fastq_stats`: Statistics for the input HiFi reads.
  - `raw_hifi_to_reference_alignment_depth`: Depth of coverage for the alignments.

- **Structural Variant Calling**

  - `raw_hifi_to_reference_alignment_structural_variants_vcf`: VCF file containing called structural variants.
  - `raw_hifi_to_reference_alignment_structural_PASS_variants_vcf`: VCF file containing structural variants that have passed all filters.
  - `raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vcf`: VCF file containing normalized structural variants.

- **Variant Calling**

  - `raw_hifi_to_reference_alignment_all_variants_vcf`: VCF file containing all called variants.
  - `raw_hifi_to_reference_alignment_all_variants_stats`: Statistics for all called variants.
  - `raw_hifi_to_reference_alignment_PASS_variants`: VCF file containing variants that have passed all filters.
  - `raw_hifi_to_reference_alignment_PASS_norm_variants`: VCF file containing normalized variants.

- **Variant Phasing**

  - `raw_hifi_to_reference_alignment_PASS_norm_phased_variants`: VCF file containing phased variants.
  - `raw_hifi_to_reference_alignment_PASS_norm_phased_stats`: Statistics for phased variants.

- **Variant Annotation**

  - `raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_variants_vcf`: VCF file containing annotated variants.
  - `raw_hifi_to_reference_alignment_PASS_norm_phased_variants_vep_stats`: Statistics for annotated variants.

- **Summary**
  - `raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary`: List of variants (tab delimited text).
  - `raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary`: List of on-target variants (tab delimited text).
  - `raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_ontarget_variants_vcf`: VCF file containing annotated on-target variants (vcf).
  - `coverage_depth_plot`: Plot of coverage depth (png).
  - `variants_summary`: Summary of all variants (tab delimited text).
  - `variants_vaf_gt0.5_summary`: Summary of all variants after filtering by variant allele frequency (VAF<0.5) (tab delimited text).
  - `sequence_summary`: Summary reads and total variants found in the sample (tab delimited text).

## Components

- **Python packages**
  - matplotlib
  - numpy
  - argparse
- **Tools**
  - pbmm2
  - seqkit
  - samtools
  - pbsv
  - bcftools
  - bedtools
  - whatshap
- **Containers**
  - ghcr.io/anand-imcm/pb-variant-call
  - google/deepvariant
  - ensemblorg/ensembl-vep
