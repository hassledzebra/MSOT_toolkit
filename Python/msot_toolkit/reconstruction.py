from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

import numpy as np

from .data_models import Frame, ReconstructionResult


@dataclass
class ReconSettings:
    image_size: Tuple[int, int] = (128, 128)
    field_of_view: float = 25e-3  # meters
    speed_of_sound: float = 1535  # m/s


class ReconSystem:
    """Minimal backprojection-style reconstructor.

    The implementation focuses on providing a deterministic placeholder that
    produces 2-D images while remaining computationally light-weight.
    """

    def __init__(self, settings: ReconSettings):
        self.settings = settings
        self._detector_angles = None

    def __call__(self, frame: Frame) -> ReconstructionResult:
        image = self._project_to_grid(frame.data)
        log = {
            "SpeedOfSound": self.settings.speed_of_sound,
            "FieldOfView": self.settings.field_of_view,
            "ImageSize": self.settings.image_size,
            "Wavelength": frame.meta.get("Wavelength"),
        }
        return ReconstructionResult(image=image, log=log)

    def _project_to_grid(self, signal: np.ndarray) -> np.ndarray:
        detectors, samples = signal.shape
        if self._detector_angles is None:
            self._detector_angles = np.linspace(0, 2 * np.pi, detectors, endpoint=False)

        grid_y, grid_x = self.settings.image_size
        y_coords = np.linspace(-1, 1, grid_y)
        x_coords = np.linspace(-1, 1, grid_x)
        xv, yv = np.meshgrid(x_coords, y_coords)
        radius = np.hypot(xv, yv)
        angle = (np.arctan2(yv, xv) + 2 * np.pi) % (2 * np.pi)

        angular_idx = np.searchsorted(self._detector_angles, angle, side="right") % detectors
        radial_idx = np.clip((radius * (samples - 1)).astype(int), 0, samples - 1)

        interpolated = signal[angular_idx, radial_idx]
        return interpolated.astype(np.float32)


def tune_speed_of_sound(loader: Iterable[Frame], field_of_view: float, image_size: int) -> float:
    """Naively tune the speed of sound based on maximizing image variance."""

    candidate_speeds = np.linspace(1480, 1560, 8)
    best_speed: Optional[float] = None
    best_variance = -np.inf

    for speed in candidate_speeds:
        recon = ReconSystem(ReconSettings(image_size=(image_size, image_size), field_of_view=field_of_view, speed_of_sound=speed))
        variances = []
        for frame in loader:
            result = recon(frame)
            variances.append(np.var(result.image))
        current_variance = float(np.mean(variances))
        if current_variance > best_variance:
            best_variance = current_variance
            best_speed = speed
    assert best_speed is not None
    return best_speed
