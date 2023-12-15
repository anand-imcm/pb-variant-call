version 1.0

# call variants using Google deep variant
task AnnotateVariants {
    
    input {
        File vcf
        File vep_cache
        File genome_reference
        String file_label
        String vep_version = "release_110.1"
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        unzip ~{vep_cache} -d vep_cache/

        perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
            --input_file ~{vcf} \
            --vcf \
            --output_file ~{file_label}raw_hifi_to_reference_alignment_PASS_norm_variants_vep_annotated.vcf \
            --stats_file ~{file_label}raw_hifi_to_reference_alignment_PASS_norm_variants_vep_stats.txt \
            --stats_text \
            --cache \
            --dir_cache vep_cache/ \
            --fasta genome_reference.fasta \
            --numbers -offline --hgvs --shift_hgvs 0 --terms SO --symbol \
            --sift b --polyphen b --total_length --ccds --canonical --biotype \
            --protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg \
            --custom file=vep_cache/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
    >>>

    output {
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_vep_annotated = file_label + "raw_hifi_to_reference_alignment_PASS_norm_variants_vep_annotated.vcf"
        File raw_hifi_to_reference_alignment_PASS_norm_phased_variants_vep_stats = file_label + "raw_hifi_to_reference_alignment_PASS_norm_variants_vep_stats.txt"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep:~{vep_version}"
        memory: "32G"
        disks: "local-disk 40 HDD"
    }
}