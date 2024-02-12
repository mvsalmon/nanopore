from argparse import ArgumentParser
from pathlib import Path

def nanoSummary(opts):
    """Produce summary file from NanoStats files"""
    # Read list of NanoStat file paths
    with open(opts.input, "r") as f:
        f_paths = f.readlines()
    
    # Main loop to extract data from each NanoStat file
    for pth in f_paths:
        pth = Path(pth.strip())
        with pth.open("r") as f:
            print(f.readlines())





if __name__ == "__main__":
    parser = ArgumentParser(description='Script to summarise NanoPlot QC results')
    parser.add_argument('-i',
                        '--input',
                        required=True,
                        #type=Path,
                        help="Path to file containing a list of paths to NanoPlot 'NanoStats' output.")
    parser.add_argument('-o',
                        '--output',
                        default='summary_output.txt',
                        help="Name of output file. Default = 'summary_output.txt'")
    opts = parser.parse_args()
    nanoSummary(opts)