"""Parsers for specie prediction tools."""
import logging

import pandas as pd

LOG = logging.getLogger(__name__)


SPP_MIN_READ_FRAC = 0.001


def parse_kraken_result(file: str):
    """parse_species_pred""Parse species prediciton result"""
    specie_pred: pd.DataFrame = pd.read_csv(file, sep="\t").sort_values(
        "fraction_total_reads", ascending=False
    )
    # limit the number of predicted species
    specie_pred = specie_pred[specie_pred["fraction_total_reads"] > SPP_MIN_READ_FRAC]
    return specie_pred.to_dict(orient="records")
