from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np


@dataclass
class Frame:
    """Single MSOT acquisition frame.

    Attributes
    ----------
    data:
        Raw time-series data. The expected shape is ``(detectors, samples)``.
    meta:
        Metadata dictionary populated by :class:`~msot_toolkit.loader.MSOTSignalLoader`.
    """

    data: np.ndarray
    meta: Dict[str, float]


@dataclass
class ReconstructionResult:
    """Result produced by :class:`~msot_toolkit.reconstruction.ReconSystem`."""

    image: np.ndarray
    log: Optional[Dict[str, str]] = None


@dataclass
class UnmixingResult:
    """Result produced by :class:`~msot_toolkit.unmixing.UnmixSystem`."""

    unmixed_image: np.ndarray
    log: Optional[Dict[str, str]] = None
