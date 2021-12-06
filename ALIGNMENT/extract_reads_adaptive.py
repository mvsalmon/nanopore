#! /usr/bin/env python

"""
extract_reads.py
Created by Tim Stuart

MS additions for handling adaptive sampling output
"""

import pysam

def get_reads(adaptive_output):
    #Added by MS
    #iterate over adaptive output file and store each read_id in dict {decision:[id1, id2, id3...]}
    read_dict = {}
    with open(adaptive_output) as f:
        #skip header row
        for head in range(1):
            next(f)
        for line in f:
            read = line.split(',')
            read_id = read[4].strip()
            decision = read[6].strip()
            if decision not in read_dict:
                read_dict[decision] = [read_id]
            else:
                read_dict[decision].append(read_id)

    return read_dict


def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n


def extract_reads(options):
    if options.names:
        n = get_names(options.names)
    elif options.adaptive_output:
        n = get_reads(options.adaptive_output)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()
    if type(n) == list:
        out = pysam.Samfile(options.out, 'wb', header=header)
        for name in n:
            try:
                name_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = name_indexed.find(name)
                for x in iterator:
                    out.write(x)
    elif type(n) == dict:
        for key in n:
            out_name = key + '_' + options.out
            out = pysam.Samfile(out_name, 'wb', header=header)
            for name in n[key]:
                try:
                    name_indexed.find(name)
                except KeyError:
                    pass
                else:
                    iterator = name_indexed.find(name)
                    for x in iterator:
                        out.write(x)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='list of read names to extract', required=False)
    parser.add_argument('-a', '--adaptive_output', help='adaptive sampling output file', required=False)
    parser.add_argument('-o', '--out', help='file name for extracted alignments', required=True)
    options = parser.parse_args()
    extract_reads(options)
