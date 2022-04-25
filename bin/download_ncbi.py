#!/usr/bin/env python3

# download_ncbi.py <REFERENCE name> to download the refseq fastq from NCBI

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

def download_ncbi(reference, reffolder):
	""" Checks available references, downloads from NCBI if not present """
	try:
	    DEVNULL = open(os.devnull, "wb")
	    Entrez.email = "2@2.com"
	    record = Entrez.efetch(
		db="nucleotide", id=reference, rettype="fasta", retmod="text"
	    )
	    sequence = record.read()
	    output = f"{reffolder}/{reference}.fasta"
	    with open(output, "w") as f:
                f.write(sequence)
	    bwaindex = f"bwa index {output}"
	    proc = subprocess.Popen(
		bwaindex.split(),
		cwd=reffolder,
		stdout=DEVNULL,
		stderr=DEVNULL,
	    )
	    out, err = proc.communicate()
	    samindex = f"samtools faidx {output}"
	    proc = subprocess.Popen(
		samindex.split(),
		cwd=reffolder,
		stdout=DEVNULL,
		stderr=DEVNULL,
	    )
	    out, err = proc.communicate()
	    print(f"Downloaded reference {reference}")
	except Exception as e:
	    print(f"Unable to download genome '{reference}' from NCBI")

parser = argparse.ArgumentParser()
parser.add_argument('reference')
parser.add_argument('download_folder')
args = parser.parse_args()
download_ncbi(args.reference, args.download_folder)
