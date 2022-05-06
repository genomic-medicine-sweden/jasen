"""Typing related data models"""

from typing import Dict, Optional

from pydantic import Field

from .base import RWModel


class ResultMlstBase(RWModel):
    """Base class for storing MLST-like typing results"""

    alleles: Dict[str, int | None]


class TypingResultMlst(ResultMlstBase):
    """MLST results"""

    scheme: str
    sequence_type: int | None = Field(None, alias="sequenceType")


class TypingResultCgMlst(ResultMlstBase):
    """MLST results"""

    n_novel: int = Field(0, alias="nNovel")
    n_missing: int = Field(0, alias="nNovel")
