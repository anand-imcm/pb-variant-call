version 1.0

# align clustered hifi reads to reference and generate the final bam using pbaa bampaint
task AlignHifiReads {
    
    input {
        File hifi_reads_fastq_gz
        File pbmm2_index
        String file_label
    }  

    String log_level = "DEBUG"

    command <<<
        set -euo pipefail

        ln -s ~{pbmm2_index} genome_reference.mmi

        gunzip -c ~{hifi_reads_fastq_gz} > ~{file_label}.hifi_reads.fastq

        pbmm2 align \
            --preset hifi \
            --sort \
            genome_reference.mmi \
            ~{file_label}.hifi_reads.fastq \
            ~{file_label}_raw_hifi_to_reference_alignment.bam \
            --rg "@RG\tID:~{file_label}\tSM:~{file_label}" \
            --log-level ~{log_level} \
            --log-file ~{file_label}_raw_hifi_to_reference_alignment.log
        
        seqkit stats -a -T ~{file_label}.hifi_reads.fastq > ~{file_label}_raw_hifi_reads_fastq_stats.tab
        
        samtools depth ~{file_label}_raw_hifi_to_reference_alignment.bam -o ~{file_label}_raw_hifi_to_reference_alignment.depth.txt

    >>>

    output {
        File raw_hifi_to_reference_alignment_bam = file_label + "_raw_hifi_to_reference_alignment.bam"
        File raw_hifi_to_reference_alignment_index = file_label + "_raw_hifi_to_reference_alignment.bam.bai"
        File raw_hifi_to_reference_alignment_log = file_label + "_raw_hifi_to_reference_alignment.log"
        File raw_hifi_reads_fastq_stats = file_label + "_raw_hifi_reads_fastq_stats.tab"
    }

}