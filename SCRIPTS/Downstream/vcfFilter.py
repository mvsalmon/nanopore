# Script to add a filter to Sniffles vcf files based on SV supporting reads

import os
from pathlib import Path
from argparse import ArgumentParser

import pandas as pd
import vcfpy

def filter_vcf(opts):
    """Adds 'MIN_SUPPORT' to FILTER field of vcf file where SUPPORT < --min_support"""
    paths = []
    # Check if input is a path to a vcf file
    sufx = Path(opts.input).suffix

    if sufx == '.vcf':
        paths.append(opts.input)
    else:
        # Read vcf file paths
        with open(opts.input, "r") as f:
            paths = f.readlines()
    
    # Set up lists to hold summary info and sample names
    summary_dfs = []
    sample_names = []
    
    # Set filter ID for header
    filter_ID = 'MIN_SUPPORT_' + str(opts.min_support)

    # Set summary file name as specified if more than 1 vcf file is being processed
    # Otherwise use name as parsed from vcf file name
    if len(paths) > 1:
        csv_output_path = os.path.join(opts.output, opts.summary_file_name)
    else:
        sample_name = os.path.basename(paths[0]).removesuffix('_sniffles.vcf')
        csv_output_path = os.path.join(opts.output, f'{sample_name}_vcf_summary.csv')
    
    # Main loop to filter input vcfs
    for vcf in paths:
        # parse file path and construct output for filtered vcf
        pth = Path(vcf.strip())
        sample_name = os.path.basename(pth).removesuffix('_sniffles.vcf')
        vcf_output_path = os.path.join(opts.output, f'{sample_name}_sniffles_filtered.vcf')  

        # Open vcf file
        reader = vcfpy.Reader.from_path(pth)

        # Add SUPPORT_MIN to filter header
        reader.header.add_filter_line(vcfpy.OrderedDict([
            ('ID', filter_ID), ('Description', f'Less than {opts.min_support} supporting reads')
            ])
            )

        # Open output file - path to output file + header object from reader     
        writer = vcfpy.Writer.from_path(vcf_output_path, reader.header)
        
        # Set up dict to hold summary info
        summary_data = {'All Variants' : {'Total': 0},
                    'Below Min. Support' : {'Total': 0},
                    }
        
        # Iterate over vcf, generate summary and store filtered output
        for record in reader:
            sv_type = record.INFO['SVTYPE']
            summary_data['All Variants']['Total'] += 1 

            # Add record to total count
            if sv_type not in summary_data['All Variants']:
                summary_data['All Variants'][sv_type] = 1
            else:
                summary_data['All Variants'][sv_type] += 1

                   
            # Count records with < --min_support
            if record.INFO['SUPPORT'] < opts.min_support:
                # Add value to FILTER field
                record.add_filter(filter_ID)
                summary_data['Below Min. Support']['Total'] += 1
                             
                if sv_type not in summary_data['Below Min. Support']:
                    summary_data['Below Min. Support'][sv_type] = 1
                else:
                    summary_data['Below Min. Support'][sv_type] += 1
            
            # Save record to filtered vcf
            writer.write_record(record)
    

        print(sample_name)
        # Make df of summary data
        sample_df = pd.DataFrame.from_dict(summary_data, orient='index')
        print(sample_df)
        summary_dfs.append(sample_df)
        sample_names.append(sample_name)
        
        
    # Join all summary dfs and save
    pd.concat(summary_dfs, keys=sample_names).to_csv(csv_output_path)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                      help='Path to file containing list of vcf files to filter, or path to single vcf file',
                      required=True)
    parser.add_argument('-o',
                        '--output',
                      help='Path to directory to write output file',
                      default='.')
    parser.add_argument('--min_support',
                        help='Minimum number of supporting reads per SV. Default = 3',
                        type=int,
                        default = 3)
    parser.add_argument('--summary_file_name',
                        help='Name to use for summary file when filtering multiple vcfs',
                        default='vcf_summary.csv')
    opts = parser.parse_args()
    print(os.getcwd())
    filter_vcf(opts)
