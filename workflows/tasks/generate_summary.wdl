version 1.0

task Summary {
    
    input {
        File vcf
        File vcfSV
        File bed
        File depth
        File raw_hifi_reads_fastq_stats
        File raw_hifi_to_reference_alignment_log
        String file_label
    }  

    command <<<
        set -euo pipefail

        # summary using all variants
        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{vcf} > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv
        modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv

        # on-target variants vcf
        bedtools intersect -header -a ~{vcf} -b ~{bed} -wa > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf

        # summary using on-target variants
        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv
        modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv

        # summary using structural variants
        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%AD:%DP:%SAC]\n" ~{vcfSV} > ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv
        modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv

        python /scripts/plot_bam_coverage.py \
            -d ~{depth} \
            -t ~{bed} \
            -p ~{file_label}
        
        perl /scripts/report.pl \
            --fastq ~{raw_hifi_reads_fastq_stats} \
            --pbmmlog ~{raw_hifi_to_reference_alignment_log} \
            --allVariants ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv \
            --onTargetVariants ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv \
            --structuralVariants ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv \
            --prefix ~{file_label}
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_ontarget_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv"
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv"
        File coverage_depth_plot = file_label + "_coverage_depth.png"
        File variants_summary = file_label + "_variants_summary.tsv"
        File sequence_summary = file_label + "_sequence_summary.tsv"
    }
    
    runtime {
        docker: "pbvarcall:latest"
        memory: "32G"
        disks: "local-disk 40 HDD"
    }
}