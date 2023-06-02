"""Generic database objects of which several other models are based on."""
from datetime import datetime
from enum import Enum

from pydantic import BaseConfig, BaseModel


class RWModel(BaseModel):
    """Base model for read/ write operations"""

    class Config(BaseConfig):
        allow_population_by_alias = True
        allow_population_by_field_name = True
        use_enum_values = True
