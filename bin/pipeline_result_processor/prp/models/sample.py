"""Data model definition of input/ output data"""
from enum import Enum
from typing import List

from pydantic import BaseModel, Field

from .base import RWModel
from .metadata import RunMetadata
from .phenotype import PhenotypeResult, PhenotypeType
from .typing import TypingResultCgMlst, TypingResultMlst, TypingMethod

SAMPLE_ID_PATTERN = r"^[a-zA-Z1-9-_]+$"


class TaxLevel(Enum):
    P = "phylum"
    C = "class"
    O = "order"
    F = "family"
    G = "genus"
    S = "specie"


class AssemblyQc(BaseModel):
    """Assembly QC metrics."""

    total_length: int
    reference_length: int
    largest_contig: int
    n_contigs: int
    n50: int
    assembly_gc: float
    reference_gc: float
    duplication_ratio: float


class SampleQc(BaseModel):
    """Collection of sample QC."""

    assembly: AssemblyQc


class SpeciesPrediction(RWModel):
    scientific_name: str = Field(..., alias="scientificName")
    tax_id: int = Field(..., alias="taxId")
    tax_level: TaxLevel = Field(..., alias="taxLevel")
    kraken_assigned_reads: int = Field(..., alias="krakenAssignedReads")
    added_reads: int = Field(..., alias="addedReads")
    fraction_total_reads: float = Field(..., alias="fractionTotalReads")


class MethodIndex(RWModel):
    type: PhenotypeType | TypingMethod
    result: PhenotypeResult | TypingResultMlst | TypingResultCgMlst


class SampleBase(RWModel):
    """Base datamodel for sample data structure"""

    sample_id: str = Field(
        ..., alias="sampleId", min_length=3, max_length=100, regex=SAMPLE_ID_PATTERN
    )
    run_metadata: RunMetadata = Field(..., alias="runMetadata")
    qc: SampleQc
    species_prediction: List[SpeciesPrediction] = Field(..., alias="speciesPrediction")


class PipelineResult(SampleBase):
    """Input format of sample object from pipeline."""

    schema_version: int = Field(..., alias="schemaVersion", gt=0)
    # optional typing
    typing_result: List[MethodIndex] = Field(..., alias="typingResult")
    # optional phenotype prediction
    phenotype_result: List[MethodIndex] = Field(..., alias="phenotypeResult")
