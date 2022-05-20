from enum import Enum
from typing import List

from .base import RWModel


class TagType(Enum):
    VIRULENCE = "virulence"
    RESISTANCE = "resistane"
    QC = "qc"


class ResistanceTag(Enum):
    VRE = "VRE"
    ESBL = "ESBL"
    MRSA = "MRSA"


class VirulenceTag(Enum):
    PVL_ALL_POS = "pos"
    PVL_LUKS_POS = "neg/pos"
    PVL_LUKF_POS = "pos/neg"
    PVL_ALL_NEG = "neg"


class TagSeverity(Enum):
    """Defined severity classes of tags"""

    INFO = "info"
    PASSED = "pass"
    WARNING = "warning"
    DANGER = "danger"


class Tag(RWModel):
    """Tag data structure."""

    type: TagType
    label: VirulenceTag | ResistanceTag
    description: str
    severity: TagSeverity


TAG_LIST = List[Tag]
