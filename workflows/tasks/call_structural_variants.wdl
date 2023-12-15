version 1.0

task CallStructuralVariants {
    
    input {
        File raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index
        File genome_reference
        String file_label
    }  

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
            ~{file_label}raw_hifi_to_reference_alignment_structural_variants.vcf
    >>>

    output {
        File raw_hifi_to_reference_alignment_structural_variants = file_label + "raw_hifi_to_reference_alignment_structural_variants.vcf"
    }
}