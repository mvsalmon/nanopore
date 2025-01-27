# subset adaptive sampling output file by deciscion, then get list of readnames for each decision

def subset_adaptive():
    """Generate dictionary of adaptive sampling decisions, and read ids for each decision"""
    read_dict = get_reads(options.adaptive_output)
    for key, values in read_dict.items():
        with open(f'{options.output_dir}/{options.run_name}_{key}_read_ids.txt', 'w') as fp:
            fp.write('\n'.join(values))


def get_reads(adaptive_output):
    """
    Iterate over adaptive sampling output file and generate a dict with lists of read ids for each
    adaptive sampling decision: 
    {adapitive_sampling_decision:[read1, read2, read3...]}
    """
    read_dict = {}
    skipped_reads = []
    with open(adaptive_output, 'r') as f:
        #skip header row
        for header in range(1):
            next(f)
        for line in f:
            # skip null lines at eof
            if line.startswith('\x00\x00'):
                continue
            read = line.split(',')
                
            try:
                read_id = read[4].strip()
                decision = read[6].strip()
                if decision not in read_dict:
                    read_dict[decision] = [read_id]
                else:
                    read_dict[decision].append(read_id)
            # Handle reads with no associated decision
            except IndexError:
                print(f'Could not find decision for read {read_id}. Skipping...')
                skipped_reads.append(read_id)

    for key in read_dict:
        print("{0} reads = {1}".format(key, len(read_dict[key])))

    # write skipped reads
    if skipped_reads:
        print("Skipped reads:")
        print(*skipped_reads, sep='\n')
        with open(f'{options.output_dir}/{options.run_name}_skipped_reads.txt', 'w') as f:
            f.write('\n'.join(str(i) for i in skipped_reads))

    return read_dict


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-a', '--adaptive_output', help='adaptive sampling output file', required=True)
    parser.add_argument('-o', '--output_dir', help='output directory', required=True)
    parser.add_argument('-n', '--run_name', help='name of sequencing run', required=True)
    options = parser.parse_args()

    subset_adaptive()
