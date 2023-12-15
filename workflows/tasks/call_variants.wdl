version 1.0

task CallVariants {
    
    input {
        File raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index
        File genome_reference
        String file_label
        String deepvariant_version = "1.6.0"
        Int deepvariant_num_shards = 12
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{raw_hifi_to_reference_alignment_bam} ~{file_label}_raw_hifi_to_reference_alignment.bam

        ln -s ~{raw_hifi_to_reference_alignment_index} ~{file_label}_raw_hifi_to_reference_alignment.bam.bai

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
        /opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --num_shards ~{deepvariant_num_shards} \
            --ref genome_reference.fasta \
            --reads ~{file_label}_raw_hifi_to_reference_alignment.bam \
            --output_vcf  ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz

        # bcftools filter PASS variants
        bcftools view -f PASS ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz
        # bcftools norm and split bialleleic sites
        bcftools norm ~{file_label}_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz -f genome_reference.fasta -m -any -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_variants.vcf.gz
    >>>

    output {
        File raw_hifi_to_reference_alignment_all_variants_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_stats = file_label + "_raw_hifi_to_reference_alignment_all_variants.visual_report.html"
        File raw_hifi_to_reference_alignment_PASS_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_PASS_norm_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_variants.vcf.gz"
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}