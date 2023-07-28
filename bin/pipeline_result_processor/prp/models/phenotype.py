"""Datamodels used for prediction results."""
from enum import Enum
from typing import Dict, List, Union

from pydantic import BaseModel, Field

from .base import RWModel


class PredictionSoftware(Enum):
    """Container for software names."""

    # phenotype
    AMRFINDER = "amrfinder"
    RESFINDER = "resfinder"
    VIRFINDER = "virulencefinder"
    MYKROBE = "mykrobe"
    SNIPPY = "snippy"
    TBPROFILER = "tbprofiler"


class VariantType(Enum):
    SUBSTITUTION = "substitution"


class ElementType(Enum):
    AMR = "AMR"
    ACID = "STRESS_ACID"
    BIOCIDE = "STRESS_BIOCIDE"
    METAL = "STRESS_METAL"
    HEAT = "STRESS_HEAT"
    VIR = "VIRULENCE"


class DatabaseReference(RWModel):
    ref_database: Union[str, None]
    ref_id: Union[str, None]


class GeneBase(BaseModel):
    """Container for gene information"""

    accession: Union[str, None]
    # prediction info
    depth: Union[float, None]
    identity: Union[float, None]
    coverage: Union[float, None]
    ref_start_pos: Union[int, None]
    ref_end_pos: Union[int, None]
    ref_gene_length: Union[int, None]
    alignment_length: Union[int, None]
    # amrfinder extra info
    contig_id: Union[str, None]
    gene_symbol: Union[str, None]
    sequence_name: Union[str, None]
    ass_start_pos: Union[int, None]
    ass_end_pos: Union[int, None]
    strand: Union[str, None]
    element_type: Union[str, None]
    element_subtype: Union[str, None]
    target_length: Union[int, None]
    res_class: Union[str, None]
    res_subclass: Union[str, None]
    method: Union[str, None]
    close_seq_name: Union[str, None]


class ResistanceGene(GeneBase, DatabaseReference):
    """Container for resistance gene information"""

    phenotypes: List[str]


class VirulenceGene(GeneBase, DatabaseReference):
    """Container for virulence gene information"""

    virulence_category: Union[str, None]


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


class ElementTypeResult(BaseModel):
    """Phenotype result data model.

    A phenotype result is a generic data structure that stores predicted genes,
    mutations and phenotyp changes.
    """

    phenotypes: Dict[str, List[str]]
    genes: List[Union[ResistanceGene, VirulenceGene]]
    mutations: List[ResistanceVariant]
