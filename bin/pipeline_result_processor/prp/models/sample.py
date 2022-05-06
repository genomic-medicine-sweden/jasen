"""Data model definition of input/ output data"""
from datetime import datetime
from email.policy import default
from enum import Enum
from multiprocessing.sharedctypes import Value
from typing import Dict, List
from unittest.mock import Base

from pydantic import BaseModel, Field, validator

from .base import DBModelMixin, ModifiedAtRWModel, RWModel
from .tags import Tag
from .typing import TypingResultCgMlst, TypingResultMlst

SAMPLE_ID_PATTERN = "^[a-zA-Z1-9-_]+$"


class TaxLevel(Enum):
    P = "phylum"
    C = "class"
    O = "order"
    F = "family"
    G = "genus"
    S = "specie"


class VariantType(Enum):
    substitution = "substitution"


class PhenotypeType(Enum):
    amr = "antimicrobial_resistance"
    chem = "chemical_resistance"
    env = "environmental_factor_resistance"
    vir = "virulence"


class SoupType(Enum):
    db = "database"
    sw = "software"


class TypingMethod(Enum):
    mlst = "mlst"
    cgmlst = "cgmlst"


class SoupVersion(RWModel):
    """Version of Software of Unknown Provenance."""

    name: str
    version: str
    type: SoupType


class RunMetadata(BaseModel):
    """Run metadata"""

    run: Dict[str, str | List[str]]
    databases: List[SoupVersion]


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


class DatabaseReference(RWModel):
    ref_database: str
    ref_id: str


class GeneBase(BaseModel):
    """Container for gene information"""

    name: str
    accession: str
    # prediction info
    depth: float | None
    identity: float
    coverage: float
    ref_start_pos: int
    ref_end_pos: int
    ref_gene_length: int
    alignment_length: int


class ResistanceGene(GeneBase, DatabaseReference):
    """Container for resistance gene information"""

    phenotypes: List[str]


class VirulenceGene(GeneBase, DatabaseReference):
    """Container for virulence gene information"""

    virulence_category: str


class VariantBase(DatabaseReference):
    """Container for mutation information"""

    variant_type: VariantType = Field(
        ..., alias="variantType"
    )  # type of mutation insertion/deletion/substitution
    genes: List[str]
    position: int
    ref_codon: str = Field(..., alias="refCodon")
    alt_codon: str = Field(..., alias="altCodon")
    # prediction info
    depth: float


class ResistanceVariant(VariantBase):
    """Container for resistance variant information"""

    phenotypes: List[str]


class PhenotypeResult(BaseModel):
    phenotypes: Dict[str, List[str]]
    genes: List[ResistanceGene | VirulenceGene]
    mutations: List[ResistanceVariant]


class SpeciesPrediction(RWModel):
    scientific_name: str = Field(..., alias="scientificName")
    tax_id: int = Field(..., alias="taxId")
    tax_level: TaxLevel = Field(..., alias="taxLevel")
    kraken_assigned_reads: int = Field(..., alias="krakenAssignedReads")
    added_reads: int = Field(..., alias="addedReads")
    fraction_total_reads: float = Field(..., alias="fractionTotalReads")


class Comment(BaseModel):
    """Contianer for comments."""

    username: str = Field(..., min_length=0)
    created_at: datetime = Field(datetime.now(), alias="createdAt")
    comment: str = Field(..., min_length=0)
    displayed: bool = True


class CommentInDatabase(Comment):
    """Comment data structure in database."""

    id: int = Field(..., alias="id")


class SampleBase(ModifiedAtRWModel):
    """Base datamodel for sample data structure"""

    sample_id: str = Field(
        ..., alias="sampleId", min_length=3, max_length=100, regex=SAMPLE_ID_PATTERN
    )
    patient_id: str | None = Field(None, alias="patientId")
    run_metadata: RunMetadata = Field(..., alias="runMetadata")
    qc: SampleQc
    species_prediction: List[SpeciesPrediction] = Field(..., alias="speciesPrediction")
    # comments and non analytic results
    comments: List[CommentInDatabase] = []
    location: str | None = Field(None, description="Location id")


class SampleInPipelineInput(SampleBase):
    """Input format of sample object from pipeline."""

    output_version: int = Field(..., alias="outputVersion", gt=0)
    # optional typing
    mlst: TypingResultMlst
    cgmlst: TypingResultCgMlst
    # optional phenotype prediction
    antimicrobial_resistance: PhenotypeResult
    chemical_resistance: PhenotypeResult
    environmental_resistance: PhenotypeResult
    virulence: PhenotypeResult


class MethodIndex(RWModel):
    type: PhenotypeType | TypingMethod
    result: PhenotypeResult | TypingResultMlst | TypingResultCgMlst


class SampleInCreate(SampleBase):
    """Basic sample information"""

    schema_version: int = Field(
        ..., alias="schemaVersion", description="Version of database schema", gt=0
    )
    # computed tags
    tags: List[Tag] = Field([])
    # optional typing
    add_typing_result: List[MethodIndex] = Field(..., alias="addTypingResult")
    # optional phenotype prediction
    add_phenotype_prediction: List[MethodIndex] = Field(
        ..., alias="addPhenotypePrediction"
    )


class SampleInDatabase(DBModelMixin, SampleInCreate):
    """Basic sample information"""

    pass
