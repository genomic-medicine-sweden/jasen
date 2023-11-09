"""Typing related data models"""

from enum import Enum
from typing import Dict, List, Optional, Union, Any

from pydantic import Field, BaseModel

from .base import RWModel


class TypingSoftware(Enum):
    """Container for software names."""

    # typing
    CHEWBBACA = "chewbbaca"
    MLST = "mlst"
    TBPROFILER = "tbprofiler"
    MYKROBE = "mykrobe"


class TypingMethod(Enum):
    MLST = "mlst"
    CGMLST = "cgmlst"
    LINEAGE = "lineage"


class LineageInformation(RWModel):
    """Base class for storing lineage information typing results"""

    lin: Union[str, None]
    family: Union[str, None]
    spoligotype: Union[str, None]
    rd: Union[str, None]
    frac: Union[str, None]
    variant: Union[str, None]
    coverage: Union[Dict, None]


class ResultMlstBase(RWModel):
    """Base class for storing MLST-like typing results"""

    alleles: Dict[str, Union[int, str, List, None]]


class ResultLineageBase(RWModel):
    """Base class for storing MLST-like typing results"""

    lineages: List[LineageInformation]#Union[int, str, None]]


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
