"""Metadata models."""
from datetime import datetime
from enum import Enum
from lib2to3.pytree import Base
from typing import Dict, List

from pydantic import BaseModel, Field

from .base import RWModel


class SoupType(Enum):
    DB = "database"
    SW = "software"


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
    workflow_name: str
    sample_name: str
    sequencing_platform: str
    sequencing_type: str
    command: str
    date: datetime


SoupVersions = List[SoupVersion]


class RunMetadata(BaseModel):
    """Run metadata"""

    run: RunInformation
    databases: SoupVersions
