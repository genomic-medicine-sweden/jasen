"""Metadata models."""
from lib2to3.pytree import Base
from pydantic import BaseModel, Field
from enum import Enum
from datetime import datetime
from typing import List, Dict
from .base import RWModel


class SoupType(Enum):
    db = "database"
    sw = "software"


class SoupVersion(BaseModel):
    """Version of Software of Unknown Provenance."""

    name: str
    version: str
    type: SoupType


class RunInformation(RWModel):
    """Store information on a run how the run was conducted."""
    pipeline: str
    version: str
    commit: str
    analysis_profile: str = Field(..., alias="analysisProfile")
    configuration_files: List[str] = Field(..., alias="configurationFiles")
    run: str
    command: str
    date: datetime


SoupVersions = List[SoupVersion]


class RunMetadata(BaseModel):
    """Run metadata"""

    run: RunInformation
    databases: SoupVersions
