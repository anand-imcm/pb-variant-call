version 1.0

import "./tasks/align_hifi.wdl" as align
import "./tasks/call_structural_variants.wdl" as svcall
import "./tasks/call_variants.wdl" as varcall
import "./tasks/annotate_variants.wdl" as annotate
import "./tasks/annotate_structural_variants.wdl" as annotateSV
import "./tasks/phase_variants.wdl" as phase
import "./tasks/generate_summary.wdl" as report

workflow main {

    String pipeline_version = "1.2.7"
    String container_src = "ghcr.io/anand-imcm/pb-variant-call:~{pipeline_version}"
    String vep_docker = "ghcr.io/anand-imcm/ensembl-vep:release_110.1"
    String deepvariant_docker = "ghcr.io/anand-imcm/deepvariant:1.5.0"

    input {
        File reads_fastq_gz
        String prefix
        File genome_ref
        File genome_index_pbmm
        File vep_cache
        File target_bed
        File region_to_plot
    }

    parameter_meta {
        reads_fastq_gz : "Input PacBio HiFi reads in .fastq.gz format."
        prefix : "Sample name. This will be used as prefix for all the output files."
        genome_ref : "Human reference genome .fasta file. Using: Homo_sapiens.GRCh38.release110.dna.chromosome.1.fa"
        genome_index_pbmm : "Reference index generated through pbmm2 in .mmi format."
        vep_cache : "VEP cache in .zip format. (Source: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). Using: Ensembl GRCh38 release v110"
        target_bed : "Coordinates for the amplified regions (target) in .bed format. (zero-based bed file without any header)"
    }

    call align.AlignHifiReads {
        input: hifi_reads_fastq_gz = reads_fastq_gz, pbmm2_index = genome_index_pbmm, file_label = prefix, docker = container_src
    }

    call varcall.CallVariants {
        input: raw_hifi_to_reference_alignment_bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam, raw_hifi_to_reference_alignment_index = AlignHifiReads.raw_hifi_to_reference_alignment_index, genome_reference = genome_ref, file_label = prefix, deepvariant_docker = deepvariant_docker
    }

    call svcall.CallStructuralVariants {
        input: raw_hifi_to_reference_alignment_bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam, raw_hifi_to_reference_alignment_index = AlignHifiReads.raw_hifi_to_reference_alignment_index, genome_reference = genome_ref, file_label = prefix, docker = container_src
    }
    
    call annotateSV.AnnotateSVs {
        input: vcf = CallStructuralVariants.raw_hifi_to_reference_alignment_structural_PASS_norm_variants, vep_cache = vep_cache, genome_reference = genome_ref, file_label = prefix, vep_docker = vep_docker
    }

    call phase.PhaseVariants {
        input: vcf = CallVariants.raw_hifi_to_reference_alignment_PASS_norm_variants, bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam, bam_index = AlignHifiReads.raw_hifi_to_reference_alignment_index, genome_reference = genome_ref, file_label = prefix, docker = container_src
    }
    
    call annotate.AnnotateVariants {
        input: vcf = PhaseVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_variants, vep_cache = vep_cache, genome_reference = genome_ref, file_label = prefix, vep_docker = vep_docker
    }

    call report.Summary {
        input: vcf = AnnotateVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_variants, vcfSV = AnnotateSVs.raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants, bed = target_bed, depth = AlignHifiReads.raw_hifi_to_reference_alignment_depth, raw_hifi_reads_fastq_stats = AlignHifiReads.raw_hifi_reads_fastq_stats, raw_hifi_to_reference_alignment_log = AlignHifiReads.raw_hifi_to_reference_alignment_log, region_to_plot = region_to_plot, file_label = prefix, docker = container_src
    }

    output {
        File raw_hifi_to_reference_alignment_bam = AlignHifiReads.raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index = AlignHifiReads.raw_hifi_to_reference_alignment_index
        File raw_hifi_to_reference_alignment_log = AlignHifiReads.raw_hifi_to_reference_alignment_log
        File raw_hifi_reads_fastq_stats = AlignHifiReads.raw_hifi_reads_fastq_stats
        File raw_hifi_to_reference_alignment_depth = AlignHifiReads.raw_hifi_to_reference_alignment_depth

        File raw_hifi_to_reference_alignment_structural_variants_vcf = CallStructuralVariants.raw_hifi_to_reference_alignment_structural_variants
        File raw_hifi_to_reference_alignment_structural_PASS_variants_vcf = CallStructuralVariants.raw_hifi_to_reference_alignment_structural_PASS_variants
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vcf = CallStructuralVariants.raw_hifi_to_reference_alignment_structural_PASS_norm_variants
        File raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants = AnnotateSVs.raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats = AnnotateSVs.raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats

        File raw_hifi_to_reference_alignment_all_variants_vcf = CallVariants.raw_hifi_to_reference_alignment_all_variants_vcf
        File raw_hifi_to_reference_alignment_all_variants_stats = CallVariants.raw_hifi_to_reference_alignment_all_variants_stats
        File raw_hifi_to_reference_alignment_PASS_variants = CallVariants.raw_hifi_to_reference_alignment_PASS_variants
        File raw_hifi_to_reference_alignment_PASS_norm_variants = CallVariants.raw_hifi_to_reference_alignment_PASS_norm_variants

        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants = PhaseVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_variants
        File raw_hifi_to_reference_alignment_PASS_norm_phased_stats = PhaseVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_stats

        File raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_variants_vcf = AnnotateVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_variants
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_vep_stats = AnnotateVariants.raw_hifi_to_reference_alignment_PASS_norm_phased_variants_vep_stats

        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary = Summary.raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary
        File raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation = Summary.raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation
        File raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_ontarget_variants_vcf = Summary.raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_ontarget_variants
        File raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary = Summary.raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary = Summary.raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary
        File raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation = Summary.raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation
        File? coverage_depth_plot = Summary.coverage_depth_plot
        File? coverage_depth_summary = Summary.coverage_depth_summary
        File variants_summary = Summary.variants_summary
        File variants_qual_gt30_summary = Summary.variants_qual_gt30_summary
        File sequence_summary = Summary.sequence_summary
    }

    meta {
        description: "A WDL-based workflow for Variant calling and annotation using PacBio HiFi reads."
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }

}