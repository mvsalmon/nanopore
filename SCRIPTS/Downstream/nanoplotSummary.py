import os
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

def nanoSummary(opts):
    """Produce summary file from NanoStats files"""
    f_paths = []

    # Check if input is a list of paths, or single file to be parsed directly.
    if opts.single_file:
        f_paths.append(opts.input)
    else:
        # Read list of NanoStat file paths
        with open(opts.input, "r") as f:
            f_paths = f.readlines()

    # Set up empty df to store run data
    run_metrics = pd.DataFrame()

    # Set output name if provided
    if opts.name:
        out_path = os.path.join(opts.output, f'{opts.name}_QC_summary.csv')

    # Main loop to extract data from each NanoStat file and merge into a single dataframe
    for pth in f_paths:
        pth = Path(pth.strip())

        # Extract sample name from file path
        s_name = os.path.basename(pth).removesuffix("NanoStats.txt")
        if not opts.name:
            out_path = os.path.join(opts.output, f'{s_name}_QC_summary.csv')
        
        # Read relevant data from file into pandas dataframe
        data = pd.read_table(pth, sep=":\s+", skiprows=1, nrows=12, names=['metric', s_name],
                             index_col=0, engine='python')
        if run_metrics.empty:
            run_metrics = data
        else:
            run_metrics = pd.merge(run_metrics, data, on='metric')

    run_metrics = run_metrics.transpose()

    run_metrics.to_csv(out_path)


if __name__ == "__main__":
    parser = ArgumentParser(description='Script to summarise NanoPlot QC results')
    parser.add_argument('-i',
                        '--input',
                        required=True,
                        #type=Path,
                        help="""Path to file containing a list of paths to NanoPlot 'NanoStats' output. Use --single_file to specify
                        input file is a single 'NanoStats' output file.""")
    parser.add_argument('-o',
                        '--output',
                        default='.',
                        help="Path of output file.")
    parser.add_argument('--single_file',
                        action="store_true",
                        help="Specify if the input file is a 'Nanostats' output.")
    parser.add_argument('-n',
                        '--name',
                        help="""Optional name for output file. If --name is not specified, output name will be generated from the file name. 
                        If a list of files is provided and --name is not specified, output name will be generated from the last file processed.""")
    opts = parser.parse_args()
    nanoSummary(opts)