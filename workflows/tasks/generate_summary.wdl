version 1.0

task Summary {
    
    input {
        File vcf
        File vcfSV
        File bed
        File region_to_plot
        File depth
        File raw_hifi_reads_fastq_stats
        File raw_hifi_to_reference_alignment_log
        String file_label
        String docker
    }  

    command <<<
        set -euo pipefail

        if [ $(grep -v "#" ~{vcf} | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation.tsv
        else
            # summary using all variants
            bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{vcf} > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv
            modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
            sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv

            headers="chrom\tpos\tref\talt\tgenotype\tgenotype_qual\tread_depth\tallele_depth\tvariant_allele_frac\tgenotype_likelihood\tVEP_Allele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE_SELECT\tMANE_PLUS_CLINICAL\tCCDS\tENSP\tRefSeq\tSOURCE\tSIFT\tPolyPhen\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_OTH_AF\tgnomADg_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tClinVar\tClinVar_CLNSIG\tClinVar_CLNREVSTAT\tClinVar_CLNDN"
            echo -e $headers > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation.tsv
            bcftools +split-vep ~{vcf} -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AD\t%VAF\t%PL]\t%CSQ\n' -d -A tab >> ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation.tsv

            # applying VAF greater than 0.5 on all variants
            bcftools view -i 'VAF>0.5' ~{vcf} > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_vep_annotated.vcf
            # summary VAF greater than 0.5 on all variants
            bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_vep_annotated.vcf > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_summary.tsv
            modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
            sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_summary.tsv

            # on-target variants vcf
            bedtools intersect -header -a ~{vcf} -b ~{bed} -wa > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf

            # summary using on-target variants
            bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv
            modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
            sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv

            # applying VAF greater than 0.5 ontarget variants
            bcftools view -i 'VAF>0.5' ~{vcf} > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_VAF_gt0.5_variants_vep_annotated.vcf
            # summary VAF greater than 0.5 ontarget variants
            bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%GQ:%DP:%AD:%VAF:%PL:%PS]\n" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_VAF_gt0.5_variants_vep_annotated.vcf > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_ontarget_variants_summary.tsv
            modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_ontarget_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
            sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_ontarget_variants_summary.tsv
        fi

        if [ $(grep -v "#" ~{vcfSV} | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv
        else
            # summary using structural variants
            bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT:%AD:%DP:%SAC]\n" ~{vcfSV} > ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv
            modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
            sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv
            sv_ann_headers="chrom\tpos\tID\tref\talt\tgenotype\tallele_depth\tread_depth\treads_support\tsv_type\tsv_end\tsv_length\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE_SELECT\tMANE_PLUS_CLINICAL\tCCDS\tENSP\tRefSeq\tSOURCE\tSIFT\tPolyPhen\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_OTH_AF\tgnomADg_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tClinVar\tClinVar_CLNSIG\tClinVar_CLNREVSTAT\tClinVar_CLNDN"
            echo -e $sv_ann_headers > ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation.tsv
            bcftools +split-vep ~{vcfSV} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t%AD\t%DP\t%SAC]\t%SVTYPE\t%END\t%SVLEN\t%CSQ\n' -d -A tab >> ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation.tsv
        fi

        if [ $(wc -l < ~{depth}) -ne 0 ]; then
            python /scripts/plot_bam_coverage.py \
                -d ~{depth} \
                -t ~{region_to_plot} \
                -p ~{file_label}
        fi
        
        perl /scripts/report.pl \
            --fastq ~{raw_hifi_reads_fastq_stats} \
            --pbmmlog ~{raw_hifi_to_reference_alignment_log} \
            --allVariants ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv \
            --onTargetVariants ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv \
            --structuralVariants ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv \
            --allVAF ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_variants_summary.tsv \
            --onTargetVAF ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_VAF_gt0.5_ontarget_variants_summary.tsv \
            --prefix ~{file_label}
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_variants_summary.tsv"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_VEP_annotation.tsv"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_annotated_ontarget_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_vep_annotated.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants_summary.tsv"
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_summary.tsv"
        File raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_VEP_annotation.tsv"
        File? coverage_depth_plot = file_label + "_coverage_depth.png"
        File variants_summary = file_label + "_variants_summary.tsv"
        File variants_vaf_gt0_5_summary = file_label + "_vaf_gt0.5_variants_summary.tsv"
        File sequence_summary = file_label + "_sequence_summary.tsv"
    }
    
    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 40 HDD"
    }
}