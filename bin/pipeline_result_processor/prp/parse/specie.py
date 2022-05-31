"""Parsers for specie prediction tools."""
import logging

import pandas as pd

from ..models.specie import SpeciePrediction, SpeciesPrediction, TaxLevel

LOG = logging.getLogger(__name__)


SPP_MIN_READ_FRAC = 0.001

TR_KRAKEN_TAX_ID = {
    "P": TaxLevel.PHYLUM,
    "C": TaxLevel.CLASS,
    "O": TaxLevel.ORDER,
    "F": TaxLevel.FAMILY,
    "G": TaxLevel.GENUS,
    "S": TaxLevel.SPECIE,
}


def parse_kraken_result(
    file: str, level: TaxLevel = TaxLevel.SPECIE
) -> SpeciesPrediction:
    """Parse species prediciton result"""
    cols = [
        "fraction_total_reads",
        "kraken_assigned_reads",
        "added_reads",
        "tax_level",
        "tax_id",
        "scientific_name",
    ]
    specie_pred: pd.DataFrame = pd.read_csv(file, sep="\t", header=None, names=cols)
    # strip indentation spaces
    specie_pred["scientific_name"] = specie_pred["scientific_name"].apply(
        lambda x: x.lstrip()
    )
    # subset prediction to desired level id
    specie_pred = specie_pred[specie_pred["tax_level"] == level.value.upper()[0]]
    # sort values on fraction total reads assigned
    specie_pred = specie_pred.sort_values("fraction_total_reads", ascending=False)
    # limit the number of predicted species
    specie_pred = specie_pred[specie_pred["fraction_total_reads"] > SPP_MIN_READ_FRAC]
    result = [
        SpeciePrediction(
            scientific_name=row["scientific_name"],
            tax_id=row["tax_id"],
            tax_level=TR_KRAKEN_TAX_ID[row["tax_level"]],
            kraken_assigned_reads=row["kraken_assigned_reads"],
            added_reads=row["added_reads"],
            fraction_total_reads=row["fraction_total_reads"],
        )
        for row in specie_pred.to_dict(orient="records")
    ]
    return result
