#!/usr/bin/env python3

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import urllib.request
import zipfile

from Bio import Entrez
import xml.etree.ElementTree as ET

description = '''
------------------------
Title: download_ncbi.py
Author(s): Markus Johansson & Ryan Kennedy
------------------------

Description:
    This script will download genomes from NCBI. Optionally it can also index and faidx the genome as well (provided the bwa and samtools executables are available).

Procedure:
    1. Search for the user-specified genome accession number.
    2. Index genome via bwa.
    3. Samtools faidx genome via samtools.

--------------------------------------------------------------------------------------------------------------------------------
'''

usage = '''
--------------------------------------------------------------------------------------------------------------------------------
Download query genomes with optional indexing.
Executed using: python3 ./download_ncbi.py [-b/--bwaidx] [-f/--faidx] -i/--input <Genome_Accesion> [<Genome_Accesion2>] -o/--output <Output_File>
--------------------------------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.0.1'
    )
parser.add_argument(
    '-b',
    '--bwaidx',
    help='run bwa idx',
    dest='bwaidx',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-f',
    '--faidx',
    help='run samtools faidx',
    dest='faidx',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-i',
    help='input file(s) (genome accession)',
    metavar='INPUT_FILE',
    dest='inputFile',
    nargs='+',
    required=True
    )
parser.add_argument(
    '-o',
    help='output folder',
    metavar='OUTPUT_FOLDER',
    dest='outputFolder',
    required=True
    )

args = parser.parse_args()

def download_ncbi(accns, reffolder, bwaidx, faidx):
	""" Checks available references, downloads from NCBI if not present """
	for accn in accns:
		try:
			DEVNULL = open(os.devnull, "wb")
			Entrez.email = "2@2.com"
			record = Entrez.efetch(db="nucleotide", id=accn, rettype="fasta", retmod="text")
			sequence = record.read()
			output = f"{reffolder}/{accn}.fasta"
			with open(output, "w+") as f:
				f.write(sequence)
		except Exception as e:
			print(f"Unable to download genome '{accn}' from NCBI")
		if bwaidx:
			try:
				bwaindex = f"bwa index {output}"
				proc = subprocess.Popen(
				bwaindex.split(),
				cwd=reffolder,
				stdout=DEVNULL,
				stderr=DEVNULL,
				)
				out, err = proc.communicate()
				print(f"bwa index of {accn} was successful")
			except FileNotFoundError:
				print(f"bwa is not installed. Please install it and try again.")
		if faidx:
			try:
				samindex = f"samtools faidx {output}"
				proc = subprocess.Popen(
				samindex.split(),
				cwd=reffolder,
				stdout=DEVNULL,
				stderr=DEVNULL,
				)
				out, err = proc.communicate()
				print(f"samtools faidx of {accn} was successful")
			except FileNotFoundError:
				print(f"samtools is not installed. Please install it and try again.")

def main():
	download_ncbi(args.inputFile, args.outputFolder, args.bwaidx, args.faidx)

if __name__ == "__main__":
	main()