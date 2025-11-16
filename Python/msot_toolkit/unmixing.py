from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .data_models import ReconstructionResult, UnmixingResult


@dataclass
class UnmixSettings:
    endmember_names: List[str] = field(default_factory=lambda: ["Hb", "HbO2"])
    spectra: Dict[str, np.ndarray] = field(default_factory=dict)

    def as_matrix(self, wavelengths: np.ndarray) -> np.ndarray:
        if not self.spectra:
            raise ValueError("No spectra provided for unmixing")
        ordered = [self.spectra[name] for name in self.endmember_names]
        matrix = np.vstack(ordered)
        if matrix.shape[1] != len(wavelengths):
            raise ValueError("Spectra length must match number of wavelengths")
        return matrix


class UnmixSystem:
    """Linear unmixing using non-negative least squares."""

    def __init__(self, settings: UnmixSettings, wavelengths: np.ndarray):
        self.settings = settings
        self.wavelengths = wavelengths
        self._endmember_matrix = self.settings.as_matrix(wavelengths)

    def __call__(self, reconstructions: List[ReconstructionResult]) -> UnmixingResult:
        stack = np.stack([frame.image for frame in reconstructions], axis=-1)
        pixels = stack.reshape(-1, stack.shape[-1])
        unmixed = self._solve_nnls(pixels)
        unmixed_image = unmixed.reshape(
            stack.shape[0], stack.shape[1], len(self.settings.endmember_names)
        )
        log = {"Endmembers": self.settings.endmember_names}
        return UnmixingResult(unmixed_image=unmixed_image.astype(np.float32), log=log)

    def _solve_nnls(self, pixels: np.ndarray) -> np.ndarray:
        import scipy.optimize

        spectra = self._endmember_matrix.T  # shape (wavelengths, endmembers)
        coefficients = np.zeros((pixels.shape[0], spectra.shape[1]), dtype=np.float32)

        for idx, spectrum in enumerate(pixels):
            result = scipy.optimize.nnls(spectra, spectrum)
            coefficients[idx] = result[0]
        return coefficients
