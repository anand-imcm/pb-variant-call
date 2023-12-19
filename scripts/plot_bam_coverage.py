import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import FuncFormatter

def parse_arguments():
    parser = argparse.ArgumentParser(prog='plot_bam_coverage.py', description='Create a coverage plot using output from samtools depth.')
    parser.add_argument('-d','--depth',type=str, help='Tab delimited file containing the depth at each position or region. (Output from: samtools depth alignment.bam > alignment.depth)', required=True)
    parser.add_argument('-t','--target',type=str, help='Target region coordinates in bed format. (zero-based bed file without any header)', required=True)
    parser.add_argument('-p','--prefix',type=str, help='Output prefix', required=True)
    return parser.parse_args()

def read_depth_file(depth_file_path):
    with open(depth_file_path, 'r') as depth_file:
        data = [line.strip().split('\t') for line in depth_file]
    return data

def read_bed_file(bed_file_path):
    bed_intervals = []
    bed_labels = []
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split('\t')
            bed_intervals.append((int(fields[1])+1, int(fields[2])))
            bed_labels.append(fields[3])
    return bed_intervals, bed_labels

def plot_coverage(positions, coverage, bed_intervals, bed_labels, plot_file, plot_title):
    smooth_positions = np.linspace(min(positions), max(positions), 300)
    spl = make_interp_spline(positions, coverage, k=3)
    smooth_coverage = spl(smooth_positions)

    plt.figure(figsize=(10, 5))
    plt.plot(smooth_positions, smooth_coverage, linestyle='-', color='b', label='Coverage')
    plt.title(plot_title)
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.grid(True)

    for i, interval in enumerate(bed_intervals):
        x_formatter = FuncFormatter(lambda x, p: format(int(x), ','))
        plt.gca().xaxis.set_major_formatter(x_formatter)
        plt.gca().add_patch(Rectangle((interval[0], 0), interval[1] - interval[0], max(smooth_coverage), color='blue', alpha=0.2))
        plt.text((interval[0] + interval[1]) / 2, max(smooth_coverage) * 0.9, bed_labels[i], ha='center', va='center', color='black')

    plt.savefig(plot_file, dpi=300)
    # plt.show()

def main():
    args = parse_arguments()
    plot_file = f'{args.prefix}_coverage_depth.png'
    plot_title = f'{args.prefix} coverage'

    data = read_depth_file(args.depth)
    positions = [int(row[1]) for row in data]
    coverage = [int(row[2]) for row in data]

    bed_intervals, bed_labels = read_bed_file(args.target)

    plot_coverage(positions, coverage, bed_intervals, bed_labels, plot_file, plot_title)

if __name__ == "__main__":
    main()