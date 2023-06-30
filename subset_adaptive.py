#subset adaptive sampling output file by deciscion, then get list of readnames for each decision

def subset_adaptive():
    read_dict = get_reads(options.adaptive_output)
    print(read_dict.keys())
    for key, values in read_dict.items():
        # out_path = options.output_dir
        # print(out_path)
        with open(f'{options.output_dir}/{key}.txt', 'w') as fp:
            fp.write('\n'.join(values))

    
def get_reads(adaptive_output):
    """iterate over adaptive output file and store each read_id in dict {decision:[id1, id2, id3...]}"""
    read_dict = {}
    read_count = 0
    with open(adaptive_output, 'r') as f:
        #skip header row
        for head in range(1):
            next(f)
        for line in f:
            read_count += 1
            read = line.split(',')
            read_id = read[4].strip()
            decision = read[6].strip()
            if decision not in read_dict:
                read_dict[decision] = [read_id]
            else:
                read_dict[decision].append(read_id)
    for n in read_dict:
        print("{0} reads = {1}".format(n, len(read_dict[n])))
    return read_dict


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-a', '--adaptive_output', help='adaptive sampling output file', required=True)
    parser.add_argument('-o', '--output_dir', help='output directory', required=True)
    options = parser.parse_args()

    subset_adaptive()
