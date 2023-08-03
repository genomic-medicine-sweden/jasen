"""Data model definition of input/ output data"""
from enum import Enum
from typing import List, Union

from pydantic import Field

from .base import RWModel
from .metadata import RunMetadata
from .phenotype import ElementTypeResult, ElementType, PredictionSoftware
from .qc import QcMethodIndex
from .typing import TypingMethod, TypingResultCgMlst, TypingResultMlst, TypingResultLineage, TypingSoftware

# disabled validation
# SAMPLE_ID_PATTERN = r"^[a-zA-Z1-9-_]+$"
# , regex=SAMPLE_ID_PATTERN


class TaxLevel(Enum):
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


class MethodIndex(RWModel):
    """Container for key-value lookup of analytical results."""

    type: Union[ElementType, TypingMethod]
    software: PredictionSoftware | TypingSoftware | None
    result: Union[ElementTypeResult, TypingResultMlst, TypingResultCgMlst, TypingResultLineage]


class SampleBase(RWModel):
    """Base datamodel for sample data structure"""

    sample_id: str = Field(..., alias="sampleId", min_length=3, max_length=100)
    run_metadata: RunMetadata = Field(..., alias="runMetadata")
    qc: List[QcMethodIndex] = Field(...)
    species_prediction: List[SpeciesPrediction] = Field(..., alias="speciesPrediction")


class PipelineResult(SampleBase):
    """Input format of sample object from pipeline."""

    schema_version: int = Field(..., alias="schemaVersion", gt=0)
    # optional typing
    typing_result: List[MethodIndex] = Field(..., alias="typingResult")
    # optional phenotype prediction
    element_type_result: List[MethodIndex] = Field(..., alias="elementTypeResult")
