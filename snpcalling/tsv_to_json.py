#!/usr/bin/env python

'''
    Usage:   python vcf_json.py <infile> <outfile>
    Example: python vcf_to_json.py transposed_report.tsv report.json
'''

import pandas
import argparse

parser = argparse.ArgumentParser(description='This script converts a \
    tab-separated report file (.tsv) from VCF \
    into json format.')

parser.add_argument('infile', help='Input file, for example \
    transposed_report.tsv')
parser.add_argument('outfile', help='Output file, chose a \
    file name, for example report.json')

args = parser.parse_args()

report = pandas.read_csv(args.infile, sep='\t', index_col=0)
#report.to_json(args.outfile, orient='index', indent=4)
 
report.reset_index(inplace=True)
report.to_json(args.outfile, orient='index', indent=4)


print('Converted {} and saved copy to {}'.format(args.infile, args.outfile))
