version 1.0

task CallStructuralVariants {
    
    input {
        File raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index
        File genome_reference
        String file_label
        String docker
        Int memory_gb = 24
        Int cpu = 16
    }

    Int disk_size_gb = ceil(size([raw_hifi_to_reference_alignment_bam, genome_reference], "GB")) * 3

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{raw_hifi_to_reference_alignment_bam} ~{file_label}_raw_hifi_to_reference_alignment.bam

        ln -s ~{raw_hifi_to_reference_alignment_index} ~{file_label}_raw_hifi_to_reference_alignment.bam.bai

        # sv call
        pbsv discover \
            ~{file_label}_raw_hifi_to_reference_alignment.bam \
            ~{file_label}_raw_hifi_to_reference_alignment.svsig.gz

        pbsv call \
            genome_reference.fasta \
            ~{file_label}_raw_hifi_to_reference_alignment.svsig.gz \
            ~{file_label}_raw_hifi_to_reference_alignment_structural_variants.vcf

        # bcftools filter PASS variants
        bcftools view -f PASS ~{file_label}_raw_hifi_to_reference_alignment_structural_variants.vcf -Ov -o ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_variants.vcf
        # bcftools norm and split bialleleic sites
        bcftools norm ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_variants.vcf -f genome_reference.fasta -m -any -Ov -o ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants.vcf
        if [ $(cat ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants.vcf| grep -v "#" | wc -l) -eq 0 ]; then
            echo "false" > annotate.txt
        else
            echo "true" > annotate.txt
        fi
    >>>

    output {
        File raw_hifi_to_reference_alignment_structural_variants = file_label + "_raw_hifi_to_reference_alignment_structural_variants.vcf"
        File raw_hifi_to_reference_alignment_structural_PASS_variants = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_variants.vcf"
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_variants.vcf"
        Boolean ann = read_boolean("annotate.txt")
    }
    
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}