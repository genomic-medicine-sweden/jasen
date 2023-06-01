"""Parsers for species prediction tools."""
import logging

import pandas as pd

LOG = logging.getLogger(__name__)


SPP_MIN_READ_FRAC = 0.001


def parse_kraken_result(file: str):
    """parse_species_pred""Parse species prediciton result"""
    tax_lvl_dict = {
        "P": "phylum",
        "C": "class",
        "O": "order",
        "F": "family",
        "G": "genus",
        "S": "species",
    }
    columns = {"name": "scientific_name"}
    species_pred: pd.DataFrame = (
        pd.read_csv(file, sep="\t")
        .sort_values("fraction_total_reads", ascending=False)
        .rename(columns=columns)
        .replace({"taxonomy_lvl": tax_lvl_dict})
    )
    # limit the number of predicted species
    species_pred = species_pred[
        species_pred["fraction_total_reads"] > SPP_MIN_READ_FRAC
    ]
    return species_pred.to_dict(orient="records")
