"""Parse output of QC tools."""

import csv
import logging
import json

from ..models.qc import QcMethodIndex, QcSoftware, QuastQcResult, PostAlignQcResult
from click.types import File

LOG = logging.getLogger(__name__)


def parse_quast_results(file: File) -> QcMethodIndex:
    """Parse quast file and extract relevant metrics.

    Args:
        sep (str): seperator

    Returns:
        AssemblyQc: list of key-value pairs
    """
    LOG.info(f"Parsing tsv file: {file.name}")
    creader = csv.reader(file, delimiter="\t")
    header = next(creader)
    raw = [dict(zip(header, row)) for row in creader]
    qc_res = QuastQcResult(
        total_length=int(raw[0]["Total length"]),
        reference_length=raw[0]["Reference length"],
        largest_contig=raw[0]["Largest contig"],
        n_contigs=raw[0]["# contigs"],
        n50=raw[0]["N50"],
        assembly_gc=raw[0]["GC (%)"],
        reference_gc=raw[0]["Reference GC (%)"],
        duplication_ratio=raw[0]["Duplication ratio"],
    )
    return QcMethodIndex(software=QcSoftware.QUAST, result=qc_res)


def parse_postalignqc_results(input_file: File) -> QcMethodIndex:
    """Parse postalignqc json file and extract relevant metrics.

    Args:
        sep (str): seperator

    Returns:
        PostAlignQc: list of key-value pairs
    """
    LOG.info(f"Parsing json file: {input_file.name}")
    qc_dict = json.load(input_file)
    qc_res = PostAlignQcResult(
        ins_size = int(float(qc_dict["ins_size"])),
        ins_size_dev = int(float(qc_dict["ins_size_dev"])),
        mean_cov = int(qc_dict["mean_cov"]),
        pct_above_x = qc_dict["pct_above_x"],
        mapped_reads = int(qc_dict["mapped_reads"]),
        tot_reads = int(qc_dict["tot_reads"]),
        iqr_median = float(qc_dict["iqr_median"]),
        dup_pct = float(qc_dict["dup_pct"]),
        dup_reads = int(qc_dict["dup_reads"]),
    )
    return QcMethodIndex(software=QcSoftware.POSTALIGNQC, result=qc_res)
