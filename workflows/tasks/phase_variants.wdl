version 1.0

task PhaseVariants {
    
    input {
        File vcf
        File bam
        File bam_index
        File genome_reference
        String file_label
        String docker
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta
        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai

        ln -s ~{bam} sample.bam
        ln -s ~{bam_index} sample.bam.bai

        # Check if the output of samtools view command has any lines
        if [ $(zcat ~{vcf}| grep -v "#" | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt
        else
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
        fi
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_stats = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt"
    }

    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 40 HDD"
    }
}