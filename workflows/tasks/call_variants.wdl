version 1.0

# call variants using Google deep variant
task CallVariantsDV {
    
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
    >>>

    output {
        File raw_hifi_to_reference_alignment_all_variants_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_stats = file_label + "_raw_hifi_to_reference_alignment_all_variants.visual_report.html"
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}