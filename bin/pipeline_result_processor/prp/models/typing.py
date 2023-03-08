"""Typing related data models"""

from enum import Enum
from typing import Dict, Optional, Union

from pydantic import Field

from .base import RWModel


class TypingMethod(Enum):
    MLST = "mlst"
    CGMLST = "cgmlst"


class ResultMlstBase(RWModel):
    """Base class for storing MLST-like typing results"""

    alleles: Dict[str, Union[int, str, None]]


class TypingResultMlst(ResultMlstBase):
    """MLST results"""

    scheme: str
    sequence_type: Union[int, None] = Field(None, alias="sequenceType")


class TypingResultCgMlst(ResultMlstBase):
    """MLST results"""

    n_novel: int = Field(0, alias="nNovel")
    n_missing: int = Field(0, alias="nNovel")
