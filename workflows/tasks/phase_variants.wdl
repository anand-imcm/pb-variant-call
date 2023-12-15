version 1.0

# call variants using Google deep variant
task PhaseVariants {
    
    input {
        File vcf
        File bam
        File bam_index
        File genome_reference
        String file_label
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta
        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai

        ln -s ~{bam} sample.bam
        ln -s ~{bam_index} sample.bam.bai

        # phasing
        whatshap phase \
            ~{vcf} \
            sample.bam \
            --reference=genome_reference.fasta \
            --indels \
            --output ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf \

        # phasing stats
        whatshap stats \
            ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf \
            --tsv=~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt

    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_stats = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt"
    }
}