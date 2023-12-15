version 1.0

task OntargetVariants {
    
    input {
        File vcf
        File bed
        String file_label
    }  

    command <<<
        set -euo pipefail

        bedtools intersect -header -a ~{vcf} -b ~{bed} -wa > ~{file_label}_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants.vcf
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants = file_label + "_raw_hifi_to_reference_alignment_PASS_norm_phased_ontarget_variants.vcf"
    }
}