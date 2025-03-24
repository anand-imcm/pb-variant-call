version 1.0

task PhaseVariants {
    
    input {
        File vcf
        File bam
        File bam_index
        File genome_reference
        String file_label
        String docker
        Int memory_gb = 24
        Int cpu = 16
    }  

    Int disk_size_gb = ceil(size([vcf, bam, genome_reference], "GB")) * 2

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta
        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai

        ln -s ~{bam} sample.bam
        ln -s ~{bam_index} sample.bam.bai

        # Check if the output of samtools view command has any lines
        if [ $(zcat ~{vcf}| grep -v "#" | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt
            echo "false" > annotate.txt
        else
            # phasing
            whatshap phase \
                ~{vcf} \
                sample.bam \
                --reference=genome_reference.fasta \
                --indels \
                --output ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf

            # phasing stats
            whatshap stats \
                ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf \
                --tsv=~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt
            echo "true" > annotate.txt
        fi
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_variants.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_stats = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_stats.txt"
        Boolean ann = read_boolean("annotate.txt")
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}