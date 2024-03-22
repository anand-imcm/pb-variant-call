# Variant Calling and Annotation using Hi-Fi Reads

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/anand-imcm/pb-variant-call/publish.yml)&nbsp;&nbsp;
![GitHub release (with filter)](https://img.shields.io/github/v/release/anand-imcm/pb-variant-call)&nbsp;&nbsp;
[![Open](https://img.shields.io/badge/Open-Dockstore-blue)](https://dockstore.org/workflows/github.com/anand-imcm/pb-variant-call:main?tab=info)

> [!TIP]
> To import the workflow into your Terra workspace, click on the above Dockstore badge, and select 'Terra' from the 'Launch with' widget on the Dockstore workflow page.

This repository contains a WDL-based workflow for variant calling and annotation using Hi-Fi reads. The workflow includes several steps such as alignment, variant calling, VCF filtering, VCF normalization, variant phasing, variant annotation, and structural variant calling.

## Workflow Steps

- **Alignment**: The HiFi reads are aligned to a reference genome using `pbmm2`. The output is a BAM file that contains the alignments.

- **Variant Calling**: Variants are called from the alignments using `deepvariant`. The output is a VCF file that contains the called variants.

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
  - `vep_cache` : VEP cache in .zip format. This [cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) is required by the `VEP` tool for annotation. The version being used is [Ensembl GRCh38 release v110](https://ftp.ensembl.org/pub/release-110/variation/vep/homo_sapiens_vep_110_GRCh38.tar.gz).
  - `target_bed` : Coordinates for the target (amplified) regions (0-based bed file). This will be used to report the on-target and off-target variants. Ensembl/Gencode gene model is currently being used.
  - `region_to_plot` : Bed file which contains all the genomic coordinates of MTX1, GBAP1, MTX1P1 and GBA1 for the depth of coverage plot.
- **optional**
  - `vep_version` : The version of the VEP tool to use. Default value: `release_110.1`. This should be compatible with the `VEP` version.
  - `deepvariant_num_shards` : The number of shards to use when running DeepVariant. Default value : `12`.
  - `deepvariant_version` : The version of the DeepVariant tool to use. Default value: `1.5.0`.

> [!NOTE]
> A custom `vep_cache.zip` file has been created which contains: [Ensembl GRCh38 release v110](https://ftp.ensembl.org/pub/release-110/variation/vep/homo_sapiens_vep_110_GRCh38.tar.gz) (extracted), [clinvar.vcf.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz) and [clinvar.vcf.gz.tbi](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi)

## Outputs

The main output files are listed below:

- **Alignment**
  - `alignment_log`: Log file for the alignment step.
  - `fastq_stats`: Statistics for the input HiFi reads.
  - `alignment_depth`: Depth of coverage for the alignments.
- **Structural Variant Calling**
  - `structural_variants_vcf`: VCF file containing called structural variants.
  - `structural_PASS_variants_vcf`: VCF file containing structural variants that have passed all filters.
  - `structural_PASS_norm_variants_vcf`: VCF file containing normalized structural variants.
- **Variant Calling**
  - `all_variants_vcf`: VCF file containing all called variants.
  - `all_variants_stats`: Statistics for all called variants.
  - `PASS_variants`: VCF file containing variants that have passed all filters.
  - `PASS_norm_variants`: VCF file containing normalized variants.
- **Variant Phasing**
  - `PASS_norm_phased_variants`: VCF file containing phased variants.
  - `PASS_norm_phased_stats`: Statistics for phased variants.
- **Variant Annotation**
  - `PASS_norm_phased_annotated_variants_vcf`: VCF file containing annotated variants.
  - `PASS_norm_phased_variants_vep_stats`: Statistics for annotated variants.
- **Summary**
  - `PASS_norm_phased_variants_summary`: List of variants (tab delimited text).
  - `PASS_norm_phased_ontarget_variants_summary`: List of on-target variants (tab delimited text).
  - `PASS_norm_phased_annotated_ontarget_variants_vcf`: VCF file containing annotated on-target variants (vcf).
  - `Structural_PASS_norm_variants_summary`: List of structural variants (tab delimited text).
  - `coverage_depth_plot`: Plot of coverage depth (png).
  - `variants_summary`: Summary of all variants (tab delimited text).
  - `variants_qual_gt30_summary`: Summary of all variants after filtering by variant QUAL > 30.
  - `sequence_summary`: Summary reads and total variants found in the sample (tab delimited text).

## Detailed Description of Annotation and Summary Outputs

### Annotation

The following annotation sources are currently being used by VEP:

- 1000 Genomes Project: phase3
- Assembly: GRCh38.p14
- COSMIC: 97
- ClinVar: 202301
- GENCODE: 44
- Genebuild: 2014-07
- HGMD-PUBLIC: 20204
- PolyPhen: 2.2.3
- Regbuild: 1.0
- SIFT: 6.2.1
- dbSNP: 154
- gnomAD exomes: r2.1.1
- gnomAD genomes: v3.1.2

The annotated VEP output is derived from the following parameters:

- `--af`: Global allele frequency from all populations in the 1000 Genomes Project to the output.
- `--af_1kg`: Allele frequency data for each of the five major populations in the 1000 Genomes Project`: African, American, East Asian, European, and South Asian.
- `--af_gnomadg`: Allele frequency data from the gnomAD genome dataset to the output.
- `--biotype`: Biotype of the transcript.
- `--canonical`: A flag for the transcripts that are marked as canonical in the Ensembl database.
- `--ccds`: Consensus CDS (CCDS) ID, if available.
- `--custom`: Custom annotation data from a ClinVar. This includes the clinical significance (CLNSIG), review status (CLNREVSTAT), and disease name (CLNDN).
- `--hgvs`: HGVS notations for the variant.
- `--mane`: A flag for the transcripts that are marked as MANE Select or MANE Plus Clinical in the Ensembl database.
- `--max_af`: Maximum allele frequency across all populations.
- `--numbers`: Exon and intron numbers.
- `--polyphen b`: PolyPhen scores, using both possible and divergent predictions.
- `--protein`: Protein sequence change.
- `--pubmed`: PubMed IDs for associated publications.
- `--shift_hgvs 0`: Disables shifting HGVS notations to account for base numbering differences between coding and genomic sequences.
- `--sift b`: SIFT scores, using both possible and divergent predictions.
- `--symbol`: Gene symbol.
- `--total_length`: Total length of the transcript.
- `--xref_refseq`: RefSeq IDs, if available.

> [!NOTE]
> The consequence annotation (`CSQ`) is added to `INFO` field of the output VCF file.
> Each variant's annotation is a string of values separated by pipe characters (|). The values include the allele, the predicted consequence, the impact, the gene symbol, the gene ID, the feature type, the feature ID, the biotype, exon and intron numbers, HGVS notations, cDNA, CDS and protein positions, amino acid changes, codon changes, existing variation IDs, distance to the feature, strand, flags, symbol source, HGNC ID, canonical and MANE flags, CCDS ID, protein ID, RefSeq ID, source, SIFT and PolyPhen scores, HGVS offset, allele frequencies, maximum allele frequency, clinical significance, somatic flag, phenotype flag, PubMed IDs, and ClinVar data.

### Summary tables

- **The sequence summary table (tsv) consists of the following information**
  - `file`: Input ID
  - `fastq_num_seqs`: Number of sequences
  - `fastq_sum_len`: Number of bases or residues, with gaps or spaces counted
  - `fastq_min_len`: Minimal sequence length, with gaps or spaces counted
  - `fastq_avg_len`: Average sequence length, with gaps or spaces counted
  - `fastq_max_len`: Maximal sequence length, with gaps or spaces counted
  - `fastq_Q1`: First quartile of sequence length, with gaps or spaces counted
  - `fastq_Q2`: Median of sequence length, with gaps or spaces counted
  - `fastq_Q3`: third quartile of sequence length, with gaps or spaces counted
  - `fastq_sum_gap`: Number of gaps
  - `fastq_N50`: [N50](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50)
  - `fastq_Q20(%)`: Percentage of bases with the quality score greater than 20
  - `fastq_Q30(%)`: percentage of bases with the quality score greater than 30
  - `fastq_GC(%)`: Percentage of GC content
  - `total_mapped_reads`: Total number of reads that were mapped to the reference genome.
  - `total_unmapped_reads`: Total number of reads that were not mapped to the reference genome.
  - `total_alignment%`: Percentage of total reads that were successfully aligned to the reference genome.
  - `total_structural_variants`: Total number of structural variants detected.
  - `total_variants`: Total number of variants (SNPs and indels) detected.
  - `total_snps`: Total number of single nucleotide polymorphisms (SNPs) detected.
  - `total_indels`: Total number of insertions and deletions (indels) detected.
  - `total_ontarget_variants`: Total number of variants detected that are within the target regions.
  - `total_ontarget_snps`: Total number of SNPs detected that are within the target regions.
  - `total_ontarget_indels`: Total number of indels detected that are within the target regions.
  - `total_variants_qual_gt30`: Total number of variants with a variant QUAL greater than 30.
  - `total_snps_qual_gt30`: Total number of SNPs with QUAL greater than 30.
  - `total_indels_qual_gt30`: Total number of indels with QUAL greater than 30.
  - `total_ontarget_variants_qual_gt30`: Total number of on-target variants with QUAL greater than 30.
  - `total_ontarget_snps_qual_gt30`: Total number of on-target SNPs with QUAL greater than 30.
  - `total_ontarget_indels_qual_gt30`: Total number of on-target indels with QUAL greater than 30.

- **The variant summary table (tsv) contains the following information**
  - `Chr`: Chromosome.
  - `Pos`: The position of the variant.
  - `Ref`: The reference base.
  - `Alt`: The alternate base.
  - `Qual`: QUAL is the Phred-scaled probability that the site has no variant. Example: `QUAL=20`: 1 % chance that there is no variant at the site. `QUAL=50`: 1 in 1e5 chance that there is no variant at the site.
  - `Filter`: : PASS if this position has passed all filters, i.e., a call is made at this position.
  - `is_on_target`: Indicates whether the variant is within the target region of interest.
  - `Sample`: Sample ID.
  - `GT:GQ:DP:AD:VAF:PL:PS`: This is taken from the FORMAT column of the VCF. Where `GT` is the Genotype, the inferred genetic state of the sample (homozygous reference, heterozygous, homozygous alternate). `GQ` is the Genotype Quality, a measure of confidence in the genotype call. `DP` is the Depth, the total number of reads covering the variant position. `AD` is the Allele Depth, the number of reads supporting each allele. `VAF` is the Variant Allele Fraction, the proportion of reads supporting the alternate allele. 
  - `PL`: Phred-scaled likelihoods for genotypes as defined in the VCF specification.
  - `PS`: Phase set, indicating variants that are in the same phased haplotype. [More details here](https://whatshap.readthedocs.io/en/latest/guide.html#phase-sets).
  
> [!NOTE]
> VAF = AD_variant / (AD_reference + AD_variant) 
> `AD`           = Allele Depth, the number of reads supporting each allele.
> `AD_reference` = Number of reads supporting the reference allele.
> `AD_variant`   = Number of reads supporting the variant allele.

> `QUAL` - quality: Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT is wrong). If ALT is ‘.’ (no variant) then this is −10log10 prob(variant), and if ALT is not ‘.’ this is −10log10 prob(no variant). If unknown, the missing value should be specified. (Numeric)

> `FILTER` - filter status: PASS if this position has passed all filters, i.e., a call is made at this position. Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. “q10;s50” might indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total number of samples. ‘0’ is reserved and should not be used as a filter String. If filters have not been applied, then this field should be set to the missing value. (String, no whitespace or semicolons permitted)

## Components

- **Python packages**
  - matplotlib
  - seaborn
  - pyarrow
  - pandas
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
