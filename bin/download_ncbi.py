#!/usr/bin/env python3

import argparse
import re
import os
import shutil
import subprocess
from zipfile import ZipFile
import requests
import time

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
    version='%(prog)s 0.1.0'
    )
parser.add_argument(
    '-c',
    '--clean',
    help='WARNING: Remove all files in download directory',
    dest='clean',
    action='store_true',
    required=False
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
    try:
        with ZipFile(zip_file, 'r') as zip_object:
            zip_object.extractall(path=outdir)
    except BadZipFile:
        print(f"⚠️  Error: Invalid or corrupt zip file {zip_file}")
        return False
    return True

def copy_file(source, destination):
    """Copy file from source to destination"""
    if os.path.exists(source):
        shutil.copy(source, destination)
        print(f"📁 Copied: {source} ➡️ {destination}")
    else:
        print(f"⚠️ Source file not found: {source}")

def mkdir(dirpath):
    """Create directory"""
    os.makedirs(dirpath, exist_ok=True)

def remove_file(filepath):
    """Remove file"""
    if os.path.exists(filepath):
        os.remove(filepath)

def remove_directory(dirpath):
    """Remove entire directory"""
    try:
        shutil.rmtree(dirpath)
    except FileNotFoundError:
        pass

def remove_dir_content(dirpath):
    """Remove everything withing directory"""
    remove_directory(dirpath)
    mkdir(dirpath)

def find_files(search_term, parent_dir):
    """Find files in given directory using regex"""
    try:
        search_files = os.listdir(parent_dir)
        return sorted([
            os.path.join(parent_dir, f) for f in search_files
            if re.search(search_term, f) and not f.endswith("~")
        ])
    except FileNotFoundError:
        print(f"⚠️  Directory does not exist: {parent_dir}")
        return []

def download_and_save_file(url, output_filepath, max_retries=5):
    """Download the file with retry and backoff"""
    headers = {
        "User-Agent": "open-source-ncbi-downloader/0.1 (+https://github.com/genomic-medicine-sweden/jasen)"
    }

    for attempt in range(max_retries):
        try:
            response = requests.get(url, stream=True, headers=headers)
            if response.status_code == 429:
                wait_time = (2 ** attempt) + 1
                print(f"🔁 Rate limited (429). Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
                continue
            response.raise_for_status()
            with open(output_filepath, 'wb') as output_file:
                for chunk in response.iter_content(chunk_size=8192):
                    output_file.write(chunk)
            print(f"✅ Downloaded: {output_filepath}")
            return True
        except requests.exceptions.RequestException as e:
            wait_time = (2 ** attempt) + 1
            print(f"⚠️  Attempt {attempt+1} failed: {e}")
            if attempt < max_retries - 1:
                print(f"⏳ Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("❌ Max retries reached. Giving up.")
    return False


def index_fasta(command_str, accn, download_dir):
    """Run indexing command"""
    try:
        proc = subprocess.Popen(
            command_str.split(),
            cwd=download_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        proc.communicate()
        print(f"🔧 Indexing complete for {accn}: {command_str}")
    except FileNotFoundError:
        print(f"❌ Indexing tool not found for {command_str}")

def download_ncbi(accns, download_dir, bwaidx, faidx):
    """Download respective NCBI files"""
    mkdir(download_dir)

    for accn in accns:
        print(f"\n🔍 Processing accession: {accn}")
        zip_path = os.path.join(download_dir, f"{accn}.zip")
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accn}/download" \
              f"?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF" \
              f"&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED&filename={accn}.zip"

        # Download the ZIP file
        if not download_and_save_file(url, zip_path):
            print(f"❌ Skipping {accn} due to download failure.")
            continue

        # Unzip and check success
        if not unzip(zip_path, download_dir):
            print(f"❌ Skipping {accn} due to unzip failure.")
            continue

        # Extract files
        try:
            data_dir = os.path.join(download_dir, f"ncbi_dataset/data/{accn}")
            fasta_files = find_files(rf'^{accn}.*\.(fna|fasta)$', data_dir)
            if not fasta_files:
                print(f"⚠️  No FASTA file found for {accn}. Skipping.")
                continue

            fasta_source = fasta_files[0]
            gff_source = os.path.join(data_dir, "genomic.gff")
            fasta_dest = os.path.join(download_dir, f"{accn}.fasta")
            gff_dest = os.path.join(download_dir, f"{accn}.gff")

            copy_file(fasta_source, fasta_dest)
            copy_file(gff_source, gff_dest)

            # Optional indexing
            if bwaidx:
                index_fasta(f"bwa index {fasta_dest}", accn, download_dir)
            if faidx:
                index_fasta(f"samtools faidx {fasta_dest}", accn, download_dir)

        except Exception as e:
            print(f"❌ Error while processing {accn}: {e}")
        finally:
            remove_directory(os.path.join(download_dir, "ncbi_dataset"))
            remove_file(zip_path)

def main():
    if args.clean:
        print(f"🧹 Cleaning output directory: {args.output_dir}")
        remove_dir_content(args.output_dir)
    download_ncbi(args.accn, args.output_dir, args.bwaidx, args.faidx)

if __name__ == "__main__":
    main()
