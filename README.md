# Variant calling and annotation using Hi-Fi reads

This WDL based workflow includes several steps such as:
- Alignment
  - The HiFi reads are aligned to a reference genome using `pbmm2`. The output is a BAM file that contains the alignments.
- Variant calling
  - Variants are called from the alignments using `DeepVariant`. The output is a VCF file that contains the called variants.
- VCF filter
  - The variants in the VCF file are filterted using `bcftools -f PASS` to include only variants that have passed all filters. 
- VCF normalization
  - The called variants are normalized using `bcftools norm`. This step ensures that all variants are represented in a standard way.
- Variant phasing
  - The "PASS" variants are phased using `whatshap phase`. The phasing is encoded in the "FORMAT" column of the VCF.
- Variant annotation
  - The "PASS" variants are annotated using `VEP` tool. The output is a VCF file that contains all the additional annotation in the "INFO" column.
- Structural variant calling
  - The structural variants are called from the alignments using `pbsv`. The output is a VCF file that contains the called structural variants.

## Input

The main inputs to the workflow are:

- `reads_fastq_gz` : "Input PacBio HiFi reads in .fastq.gz format."
- `prefix` : "Sample name. This will be used as prefix for all the output files."
- `genome_ref` : "Human reference genome .fasta file."
- `genome_index_pbmm` : "Reference index generated through pbmm2 in .mmi format."
- `vep_cache` : "VEP cache in .zip format. (Source: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)"

## Output

The main outputs from the workflow are:

- A BAM file that contains the alignments of the HiFi reads to the reference genome.
- A report from `qualimap bamqc` that contains quality control metrics for the BAM file.
- A VCF file that contains the called variants.
- A VCF file that contains the normalized variants.
- A VCF file that contains the called structural variants.
- A VCF file that contains the phased variants.
- A text file that contains statistics about the phasing.

These files are placed in the `output_dir` directory.
