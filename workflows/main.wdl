version 1.0

import "./tasks/align_hifi.wdl" as align
import "./tasks/call_variants.wdl" as varcall
import "./tasks/annotate_variants.wdl" as annotate

workflow main {

    String pipeline_version = "1.0.0"

    input {
        File reads_fastq_gz
        String prefix
        File genome_ref
        File genome_index_pbmm
        File vep_cache
    }

    parameter_meta {
        reads_fastq_gz : "Input PacBio HiFi reads in .fastq.gz format."
        prefix : "Sample name. This will be used as prefix for all the output files."
        genome_ref : "Human reference genome .fasta file."
        genome_index_pbmm : "Reference index generated through pbmm2 in .mmi format."
        vep_cache : "VEP cache in .zip format. (Source: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)"
    }

    call align.AlignHifiReads {
        input: hifi_reads_fastq_gz = reads_fastq_gz, pbmm2_index = genome_index_pbmm, file_label = prefix
    }

    call varcall.CallVariantsDV {
        input: raw_hifi_to_reference_alignment_bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam, raw_hifi_to_reference_alignment_index = AlignHifiReads.raw_hifi_to_reference_alignment_index, genome_reference = genome_ref, file_label = prefix
    }
    
    call annotate.AnnotateVariantsVEP {
        input: variants_vcf = CallVariantsDV.raw_hifi_to_reference_alignment_PASS_norm_variants, vep_cache = vep_cache, genome_reference = genome_ref, file_label = prefix
    }

    output {
        File raw_hifi_to_reference_alignment_bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index = AlignHifiReads.raw_hifi_to_reference_alignment_index
        File raw_hifi_to_reference_alignment_log = AlignHifiReads.raw_hifi_to_reference_alignment_log
        File raw_hifi_reads_fastq_stats = AlignHifiReads.raw_hifi_reads_fastq_stats

        File raw_hifi_to_reference_alignment_all_variants_vcf = CallVariantsDV.raw_hifi_to_reference_alignment_all_variants_vcf
        File raw_hifi_to_reference_alignment_all_variants_stats = CallVariantsDV.raw_hifi_to_reference_alignment_all_variants_stats
        File raw_hifi_to_reference_alignment_PASS_variants = CallVariantsDV.raw_hifi_to_reference_alignment_PASS_variants
        File raw_hifi_to_reference_alignment_PASS_norm_variants = CallVariantsDV.raw_hifi_to_reference_alignment_PASS_norm_variants

        File raw_hifi_to_reference_alignment_all_variants_vep = AnnotateVariantsVEP.raw_hifi_to_reference_alignment_all_variants_vep
        File raw_hifi_to_reference_alignment_all_variants_vep_stats = AnnotateVariantsVEP.raw_hifi_to_reference_alignment_all_variants_vep_stats
    }

    meta {
        description: "A WDL-based workflow for Variant calling and annotation using PacBio HiFi reads."
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }

}