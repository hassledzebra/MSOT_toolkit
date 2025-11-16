from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.signal import butter, filtfilt

from .data_models import Frame


@dataclass
class PreFilterSettings:
    low_cut: float = 100e3
    high_cut: float = 8e6
    sampling_rate: float = 40e6
    order: int = 4


class MSOTPreFilter:
    """Simple bandpass filter to mimic the MATLAB Wiener prefilter."""

    def __init__(self, settings: PreFilterSettings):
        self.settings = settings
        self._b, self._a = self._design_filter()

    def _design_filter(self) -> Tuple[np.ndarray, np.ndarray]:
        nyquist = 0.5 * self.settings.sampling_rate
        normalized = [self.settings.low_cut / nyquist, self.settings.high_cut / nyquist]
        return butter(self.settings.order, normalized, btype="bandpass")

    def __call__(self, frame: Frame) -> Frame:
        filtered = filtfilt(self._b, self._a, frame.data, axis=-1)
        return Frame(filtered.astype(np.float32), frame.meta)
