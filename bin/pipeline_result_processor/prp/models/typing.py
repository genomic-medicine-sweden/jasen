"""Typing related data models"""

from enum import Enum
from typing import Dict, Optional, Union, Any

from pydantic import Field, BaseModel

from .base import RWModel


class TypingSoftware(Enum):
    """Container for software names."""

    # typing
    CHEWBBACA = "chewbbaca"
    MLST = "mlst"
    TBPROFILER = "tbprofiler"


class TypingMethod(Enum):
    MLST = "mlst"
    CGMLST = "cgmlst"
    LINEAGE = "lineage"


class ResultMlstBase(RWModel):
    """Base class for storing MLST-like typing results"""

    alleles: Dict[str, Union[int, str, None]]

class ResultLineageBase(RWModel):
    """Base class for storing MLST-like typing results"""

    lineages: Dict[str, Any]#Union[int, str, None]]


class TypingResultMlst(ResultMlstBase):
    """MLST results"""

    scheme: str
    sequence_type: Union[int, None] = Field(None, alias="sequenceType")


class TypingResultCgMlst(ResultMlstBase):
    """MLST results"""

    n_novel: int = Field(0, alias="nNovel")
    n_missing: int = Field(0, alias="nNovel")

class TypingResultLineage(ResultLineageBase):
    """Lineage results"""

    main_lin: str
    sublin: str