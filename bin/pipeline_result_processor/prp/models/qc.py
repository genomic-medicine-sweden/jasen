"""QC data models."""
from enum import Enum

from pydantic import BaseModel

from .base import RWModel


class QcTool(Enum):
    """Valid tools."""

    QUAST = "quast"
    FASTQC = "fastqc"


class QuastQcResult(BaseModel):
    """Assembly QC metrics."""

    total_length: int
    reference_length: int
    largest_contig: int
    n_contigs: int
    n50: int
    assembly_gc: float
    reference_gc: float
    duplication_ratio: float


class QcMethodIndex(RWModel):
    tool: QcTool
    version: str | None
    result: QuastQcResult
