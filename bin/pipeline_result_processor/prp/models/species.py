from enum import Enum
from typing import List

from pydantic import Field

from .base import RWModel


class TaxLevel(Enum):
    """Braken phylogenetic level."""

    P = "phylum"
    C = "class"
    O = "order"
    F = "family"
    G = "genus"
    S = "species"


class SpeciesPrediction(RWModel):
    scientific_name: str = Field(..., alias="scientificName")
    taxonomy_id: int = Field(..., alias="taxId")
    taxonomy_lvl: TaxLevel = Field(..., alias="taxLevel")
    kraken_assigned_reads: int = Field(..., alias="krakenAssignedReads")
    added_reads: int = Field(..., alias="addedReads")
    fraction_total_reads: float = Field(..., alias="fractionTotalReads")


SpeciesPrediction = List[SpeciesPrediction]
