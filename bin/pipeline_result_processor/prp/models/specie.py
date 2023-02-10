from enum import Enum
from typing import List

from pydantic import Field

from .base import RWModel


class TaxLevel(Enum):
    """Braken phylogenetic level."""

    PHYLUM = "P"
    CLASS = "C"
    ORDER = "O"
    FAMILY = "F"
    GENUS = "G"
    SPECIE = "S"


class SpeciePrediction(RWModel):
    name: str = Field(..., alias="scientificName")
    taxonomy_id: int = Field(..., alias="taxId")
    taxonomy_lvl: TaxLevel = Field(..., alias="taxLevel")
    kraken_assigned_reads: int = Field(..., alias="krakenAssignedReads")
    added_reads: int = Field(..., alias="addedReads")
    fraction_total_reads: float = Field(..., alias="fractionTotalReads")


SpeciesPrediction = List[SpeciePrediction]
