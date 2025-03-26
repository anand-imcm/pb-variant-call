version 1.0

task CallVariants {
    
    input {
        File raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index
        File genome_reference
        String file_label
        String deepvariant_docker = "google/deepvariant:1.5.0"
        Int deepvariant_num_shards = 12
        Int memory_gb = 24
        Int cpu = 16
    }  

    Int disk_size_gb = ceil(size([raw_hifi_to_reference_alignment_bam, genome_reference], "GB")) * 3

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{raw_hifi_to_reference_alignment_bam} ~{file_label}_raw_hifi_to_reference_alignment.bam

        ln -s ~{raw_hifi_to_reference_alignment_index} ~{file_label}_raw_hifi_to_reference_alignment.bam.bai

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai

        # Check if the output of samtools view command has any lines
        if [ $(samtools view ~{file_label}_raw_hifi_to_reference_alignment.bam | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz ~{file_label}_raw_hifi_to_reference_alignment_all_variants.visual_report.html ~{file_label}_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_variants.vcf.gz
        else
            /opt/deepvariant/bin/run_deepvariant \
                --model_type PACBIO \
                --vcf_stats_report=true \
                --num_shards $(nproc) \
                --ref genome_reference.fasta \
                --reads ~{file_label}_raw_hifi_to_reference_alignment.bam \
                --output_vcf ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz

            # bcftools filter PASS variants
            bcftools view -f PASS ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz
            # bcftools norm and split bialleleic sites
            bcftools norm ~{file_label}_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz -f genome_reference.fasta -m -any -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_variants.vcf.gz
        fi
    >>>

    output {
        File raw_hifi_to_reference_alignment_all_variants_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_stats = file_label + "_raw_hifi_to_reference_alignment_all_variants.visual_report.html"
        File raw_hifi_to_reference_alignment_PASS_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_PASS_norm_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_variants.vcf.gz"
    }

    runtime {
        docker: "~{deepvariant_docker}"
        bootDiskSizeGb: 15
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}