from enum import Enum
from typing import List

from .base import RWModel


class TagType(Enum):
    virulence = "virulence"
    resistance = "resistane"
    qc = "qc"


class ResistanceTag(Enum):
    vre = "VRE"
    esbl = "ESBL"
    mrsa = "MRSA"


class VirulenceTag(Enum):
    pvl_all_pos = "pos"
    pvl_lukS_pos = "neg/pos"
    pvl_lukF_pos = "pos/neg"
    pvl_all_neg = "neg"


class TagSeverity(Enum):
    """Defined severity classes of tags"""

    info = "info"
    passed = "pass"
    warning = "warning"
    danger = "danger"


class Tag(RWModel):
    """Tag data structure."""

    type: TagType
    label: VirulenceTag | ResistanceTag
    description: str
    severity: TagSeverity


TAG_LIST = List[Tag]
