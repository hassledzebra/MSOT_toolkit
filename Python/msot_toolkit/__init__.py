"""Lightweight Python port of the MATLAB MSOT toolkit."""

from .data_models import Frame, ReconstructionResult, UnmixingResult
from .loader import MSOTSignalLoader, create_synthetic_dataset
from .filtering import MSOTPreFilter
from .reconstruction import ReconSystem
from .unmixing import UnmixSystem
from .pipeline import run_pipeline, PipelineSettings

__all__ = [
    "Frame",
    "ReconstructionResult",
    "UnmixingResult",
    "MSOTSignalLoader",
    "create_synthetic_dataset",
    "MSOTPreFilter",
    "ReconSystem",
    "UnmixSystem",
    "PipelineSettings",
    "run_pipeline",
]
