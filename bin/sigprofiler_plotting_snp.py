#!/usr/bin/env python3

from sigProfilerPlotting import plotSBS, plotDBS
from pathlib import Path

import argparse, sys


def plot_spectra(matrix_dir, project, output_dir):
    matrix_dir = Path(matrix_dir)
    # SBS 
    matrix_dir_sbs = matrix_dir / "SBS"
    matrix_paths = matrix_dir_sbs.glob("*.all")
    output_path = f'{output_dir}/SBS/'
    for matrix_path in matrix_paths:
        plot_type=matrix_path.name.split('.')[1][3:]
        plotSBS(str(matrix_path), output_path, project, plot_type, percentage=False)
        
    # DBS
    matrix_dir_dbs = matrix_dir / "DBS"
    matrix_paths = matrix_dir_dbs.glob("*.all")
    output_path = f'{output_dir}/DBS/'
    for matrix_path in matrix_paths:
        plot_type=matrix_path.name.split('.')[1][3:]
        plotDBS(str(matrix_path), output_path, project, plot_type, percentage=False)

    
def get_arguments():
    parser = argparse.ArgumentParser(description='plot SBS and DBS spectra for SNPs')
    parser.add_argument('--matrix_dir', required=True,
                        help='path to the directory that contains the matrices, ie output from SigProfilerMatrixGenerator', type = str)
    parser.add_argument('--project', required=True,
                        help='path to the concordance file', type = str)
    parser.add_argument('--output_path', '-o', required=True,
                        help='output directory, has to have a slash for sigprofiler plotting to work', type = str)
    return parser

def main(args):
    plot_spectra(args.matrix_dir, args.project, args.output_path)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))

            
