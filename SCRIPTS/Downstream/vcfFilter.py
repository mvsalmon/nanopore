# Script to add a filter to Sniffles vcf files base on SV supporting reads
import os
from pathlib import Path
from argparse import ArgumentParser

import pandas as pd
import vcfpy

def filter_vcf(opts):
    """Adds 'SUPPORT_MIN' to FILTER field of vcf file where SUPPORT < --min_support"""
    # Read vcf file paths
    with open(opts.input, "r") as f:
        paths = f.readlines()
    
    # Set up lists to hold summary info and sample names
    summary_dfs = []
    sample_names = []
    
    for vcf in paths:
        # parse file path and construct output for filtered vcf and summary csv
        pth = Path(vcf.strip())
        sample_name = os.path.basename(pth).removesuffix('_sniffles.vcf')
        vcf_output_path = os.path.join(opts.output, f'{sample_name}_sniffles_filtered.vcf')
        csv_output_path = os.path.join(opts.output, f'{sample_name}_vcf_summary.csv')

        # Open vcf file
        reader = vcfpy.Reader.from_path(pth)

        # Open output file - path to output file + header object from reader
        writer = vcfpy.Writer.from_path(vcf_output_path, reader.header)
        
        # Set up dict to hold summary info
        summary_data = {'Total' : {'Vars': 0},
                    'Min. Support' : {'Vars': 0},
                    }
        
        # Iterate over vcf, generate summary and store filtered output
        for record in reader:
            sv_type = record.INFO['SVTYPE']
            summary_data['Total']['Vars'] += 1 
           
            # Records with < --min_support
            if record.INFO['SUPPORT'] < opts.min_support:
                # Add value to FILTER field
                record.add_filter(f'SUPPORT_MIN_{opts.min_support}')
                             
                if sv_type not in summary_data['Total']:
                    summary_data['Total'][sv_type] = 1
                else:
                    summary_data['Total'][sv_type] += 1
            # Records with > --min_support
            else:
                summary_data['Min. Support']['Vars'] += 1
                if sv_type not in summary_data['Min. Support']:
                    summary_data['Min. Support'][sv_type] = 1
                else:
                    summary_data['Min. Support'][sv_type] += 1 
            
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
                      help='Path to file containing list of vcf files to filter',
                      required=True)
    parser.add_argument('-o',
                        '--output',
                      help='Path to directory to write output file',
                      default='.')
    parser.add_argument('--min_support',
                        help='Minimum number of supporting reads per SV',
                        type=int,
                        required=True)
    opts = parser.parse_args()
    print(os.getcwd())
    filter_vcf(opts)
