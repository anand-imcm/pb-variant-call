version 1.0

task AnnotateSVs {
    
    input {
        File vcf
        File vep_cache
        File genome_reference
        String file_label
        String vep_docker = "ensemblorg/ensembl-vep:release_110.1"
        Int vep_fork = 12
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        # Check if the output of samtools view command has any lines
        if [ $(grep -v "#" ~{vcf} | wc -l) -eq 0 ]; then
            touch ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants.vcf ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats.txt
        else
            unzip ~{vep_cache} -d vep_cache/

            perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
                --input_file ~{vcf} \
                --vcf \
                --output_file ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants.vcf \
                --stats_file ~{file_label}_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats.txt \
                --stats_text \
                --cache \
                --dir_cache vep_cache/ \
                --fasta genome_reference.fasta \
                --fork ~{vep_fork} \
                --numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol \
                --sift b --polyphen b --total_length --ccds --canonical --biotype \
                --protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg \
                --custom file=vep_cache/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
        fi
    >>>

    output {
        File raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_annotated_variants.vcf"
        File raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats = file_label + "_raw_hifi_to_reference_alignment_structural_PASS_norm_variants_vep_stats.txt"
    }

    runtime {
        docker: "~{vep_docker}"
        memory: "32G"
        disks: "local-disk 40 HDD"
    }
}