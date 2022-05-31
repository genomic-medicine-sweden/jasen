from enum import Enum
from typing import List

from pydantic import Field

from .base import RWModel


class TaxLevel(Enum):
    """Braken phylogenetic level."""

    PHYLUM = "phylum"
    CLASS = "class"
    ORDER = "order"
    FAMILY = "family"
    GENUS = "genus"
    SPECIE = "specie"


class SpeciePrediction(RWModel):
    scientific_name: str = Field(..., alias="scientificName")
    tax_id: int = Field(..., alias="taxId")
    tax_level: TaxLevel = Field(..., alias="taxLevel")
    kraken_assigned_reads: int = Field(..., alias="krakenAssignedReads")
    added_reads: int = Field(..., alias="addedReads")
    fraction_total_reads: float = Field(..., alias="fractionTotalReads")


SpeciesPrediction = List[SpeciePrediction]
