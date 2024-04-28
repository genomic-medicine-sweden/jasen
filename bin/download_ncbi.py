#!/usr/bin/env python3

import argparse
import re
import os
import shutil
import subprocess
from zipfile import ZipFile
import requests

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
    '-i', '--accn',
    help='input accession number(s) (genome accession)',
    metavar='ACCN_NUMBER',
    dest='accn',
    nargs='+',
    required=True
    )
parser.add_argument(
    '-o',
    help='output directory',
    metavar='OUTPUT_DIRECTORY',
    dest='output_dir',
    required=True
    )

args = parser.parse_args()

def unzip(zip_file, outdir):
    """Unzip zip file"""
    with ZipFile(zip_file, 'r') as zip_object:
        zip_object.extractall(path=outdir)

def copy_file(source, destination):
    """Copy file from source to destination"""
    try:
        shutil.copy(source, destination)
        print(f"File copied from {source} to {destination}")
    except Exception as error_code:
        print(f"Error copying file: {error_code}")

def remove_file(filepath):
    """Remove file."""
    if os.path.exists(filepath):
        os.remove(filepath)
        print(f"The file '{filepath}' was successfully removed.")

def remove_directory(directory):
    """Remove an entire directory."""
    try:
        shutil.rmtree(directory)
        print(f"Directory '{directory}' successfully removed.")
    except FileNotFoundError:
        print(f"Directory '{directory}' does not exist.")
    except PermissionError:
        print(f"Permission denied to remove directory '{directory}'.")
    except Exception as e:
        print(f"Error occurred while removing directory '{directory}': {e}")

def find_files(search_term, parent_dir):
    """Find files in given directory using regex search term"""
    try:
        search_files = os.listdir(parent_dir)
    except FileNotFoundError:
        print(f"WARN: {parent_dir} does not exist! Trying to fix.")
    finally:
        search_files = os.listdir(parent_dir)
        found_files = sorted([os.path.join(parent_dir, search_file)
                                for search_file in search_files
                                if re.search(search_term, search_file) and
                                not search_file.endswith("~")
                            ])
        return found_files

def download_and_save_file(url, output_filepath):
    """Download the file and save it to the user-specified path"""
    try:
        # Make a request to the URL
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an error for bad responses

        # Open the output file in binary write mode
        with open(output_filepath, 'wb') as output_file:
            # Iterate over the content in chunks and write to the file
            for chunk in response.iter_content(chunk_size=8192):
                output_file.write(chunk)

        print(f"File downloaded and saved to: {output_filepath}")

    except requests.exceptions.RequestException as error_code:
        print(f"Error downloading the file: {error_code}")

def index_fasta(accn, idx_cmd, download_dir, DEVNULL):
    """Index fasta file"""
    try:
        bwaindex = idx_cmd
        proc = subprocess.Popen(
        bwaindex.split(),
        cwd=download_dir,
        stdout=DEVNULL,
        stderr=DEVNULL,
        )
        out, err = proc.communicate()
        print(f"Indexing of {accn} was successful")
    except FileNotFoundError:
        print(f"ERROR: Software to perform indexing is not installed. Please install it and try again.")

def download_ncbi(accns, download_dir, bwaidx, faidx):
    """Download gff of genome genes"""
    DEVNULL = open(os.devnull, "wb")
    for accn in accns:
        output_fpath = os.path.join(download_dir, f"{accn}.fasta")
        if bwaidx:
            index_fasta(f"bwa index {output_fpath}", accn, download_dir, DEVNULL)
        if faidx:
            index_fasta(f"samtools faidx {output_fpath}", accn, download_dir, DEVNULL)
        zip_filepath = os.path.join(download_dir, f"{accn}.zip")
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accn}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED&filename={accn}.zip"
        try:
            download_and_save_file(url, zip_filepath)
            unzip(zip_filepath, download_dir)
            gff_sourcepath = os.path.join(download_dir, f"ncbi_dataset/data/{accn}/genomic.gff")
            fasta_sourcedir = os.path.join(download_dir, f"ncbi_dataset/data/{accn}")
            fasta_sourcepath = find_files(rf'^{accn}.*\.(fna|fasta)$', fasta_sourcedir)[0]
            gff_destpath = os.path.join(download_dir, f"{accn}.gff")
            # fasta_destpath = os.path.join(download_dir, os.path.basename(fasta_sourcepath))
            fasta_destpath = os.path.join(download_dir, f"{accn}.fasta")
            copy_file(gff_sourcepath, gff_destpath)
            copy_file(fasta_sourcepath, fasta_destpath)
            remove_directory(os.path.join(download_dir, "ncbi_dataset"))
            remove_file(os.path.join(download_dir, f"{accn}.zip"))
            
        except Exception as error_code:
            print(f"Error downloading the gff file: {error_code}")
            
def download_ncbi_fasta(accns, download_dir, bwaidx, faidx, db="nucleotide"):
    """ Checks available references, downloads from NCBI if not present """
    for accn in accns:
        try:
            DEVNULL = open(os.devnull, "wb")
            Entrez.email = "2@2.com"
            record = Entrez.efetch(db=db, id=accn, rettype="fasta", retmod="text")
            sequence = record.read()
            output = f"{download_dir}/{accn}.fasta"
            with open(output, "w+") as f:
                f.write(sequence)
        except Exception as e:
            print(f"Unable to download genome '{accn}' from NCBI")
        if bwaidx:
            try:
                bwaindex = f"bwa index {output}"
                proc = subprocess.Popen(
                bwaindex.split(),
                cwd=download_dir,
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
                cwd=download_dir,
                stdout=DEVNULL,
                stderr=DEVNULL,
                )
                out, err = proc.communicate()
                print(f"samtools faidx of {accn} was successful")
            except FileNotFoundError:
                print(f"samtools is not installed. Please install it and try again.")

def main():
    if any(accn.startswith("NC") for accn in args.accn):
        download_ncbi_fasta(args.accn, args.output_dir, args.bwaidx, args.faidx)
    else:
        download_ncbi(args.accn, args.output_dir, args.bwaidx, args.faidx)

if __name__ == "__main__":
    main()
