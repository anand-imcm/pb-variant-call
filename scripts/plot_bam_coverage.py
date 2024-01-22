import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from matplotlib.ticker import FuncFormatter
import seaborn as sns

def parse_arguments():
    parser = argparse.ArgumentParser(prog='plot_bam_coverage.py', description='Create a coverage plot using output from samtools depth.')
    parser.add_argument('-d','--depth',type=str, help='Tab delimited file containing the depth at each position or region. (Output from: samtools depth alignment.bam > alignment.depth)', required=True)
    parser.add_argument('-t','--target',type=str, help='Target region coordinates in bed format. (zero-based bed file without any header) This bed file will be used to highlight the genomic region on the plot.', required=True)
    parser.add_argument('-p','--prefix',type=str, help='Output prefix', required=True)
    return parser.parse_args()

def draw_genes_and_coverage(targets_df, cov_df, plot_file, plot_title):
    fig, (ax_genes, ax_coverage) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [1, 3]}, sharex=True)
    targets_df = targets_df.sort_values(by="Start")
    chromosome_positions = {}

    for index, row in targets_df.iterrows():
        chrom, start, end, gene_info, _, strand = row["Chr"], row["Start"], row["End"], row["Label"], row["Info"], row["Strand"]
        
        start, end = int(start), int(end)

        if chrom not in chromosome_positions:
            chromosome_positions[chrom] = 0
        y_position = chromosome_positions[chrom]

        face_color = '#ADD8E6' if gene_info.startswith('GBA') else '#FFC0CB'

        rect = patches.Rectangle((start, y_position), end - start, 1, linewidth=1, edgecolor='black', facecolor=face_color)
        
        ax_genes.add_patch(rect)

        ax_genes.text(start + (end - start) / 2, y_position + 0.5, gene_info, ha='center', va='center', fontsize=8, color='black')

        chromosome_positions[chrom] += 1

    ax_genes.set_title("Region")
    x_axis_min = min(targets_df["Start"].min(), cov_df["Pos"].min()) - 5000
    x_axis_max = max(targets_df["End"].max(), cov_df["Pos"].max()) + 5000
    ax_genes.set_xlim(x_axis_min, x_axis_max)
    ax_genes.set_ylim(0, max(chromosome_positions.values()) + 0.5)  # Adjust y-axis limits as needed
    ax_genes.set_ylabel('Genes')
    ax_genes.yaxis.set_visible(False)

    ax_coverage.fill_between(cov_df['Pos'], cov_df['Depth'], alpha=0.3, color='gray')
    ax_coverage.set_title(plot_title)
    ax_coverage.set_xlabel("Position")
    ax_coverage.set_ylabel("Depth")
    ax_coverage.set_xlim(x_axis_min, x_axis_max)
    ax_coverage.set_ylim(0, max(cov_df['Depth']))
    
    x_formatter = FuncFormatter(lambda x, p: format(int(x), ','))
    ax_coverage.xaxis.set_major_formatter(x_formatter)

    plt.savefig(plot_file, dpi=300)

def main():
    args = parse_arguments()
    
    plot_file = f'{args.prefix}_coverage_depth.png'
    
    plot_title = f'{args.prefix} coverage'

    targets = pd.read_csv(args.target, sep='\t', header=None, names=["Chr", "Start", "End", "Label", "Info", "Strand"])
    
    coverage = pd.read_csv(args.depth, sep='\t', header=None, names=["Chr", "Pos", "Depth"])

    draw_genes_and_coverage(targets, coverage, plot_file, plot_title)

if __name__ == "__main__":
    main()