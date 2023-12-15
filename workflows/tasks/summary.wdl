version 1.0

task GenerateSummary {
    
    input {
        File depth
        File bed
        String file_label
    }  

    command <<<
        set -euo pipefail

        python /home/anand/Documents/aspire-files/data-oxford/terra.bio/pb-variant-call/scripts/plot_bam_coverage.py \
            -d ~{depth} \
            -t ~{bed} \
            -p ~{file_label}
    >>>

    output {
        File coverage_depth_plot = file_label + "_coverage_depth.png"
    }
}