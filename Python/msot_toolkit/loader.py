from __future__ import annotations

import pathlib
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List

import h5py
import numpy as np

from .data_models import Frame


@dataclass
class LoaderMeta:
    """Metadata describing the dataset and its acquisition."""

    wavelengths: np.ndarray
    repetition: int
    detector_count: int
    sampling_rate: float

    @property
    def frame_count(self) -> int:
        return len(self.wavelengths) * self.repetition


class MSOTSignalLoader(Iterable[Frame]):
    """Load MSOT signal data from ``.msot``/HDF5/``.npz`` files.

    The implementation is intentionally minimal. It focuses on providing a
    consistent API that mirrors the MATLAB loader while remaining easy to mock
    during development.
    """

    def __init__(self, file_path: str | pathlib.Path):
        self.file_path = pathlib.Path(file_path)
        self._frames: List[Frame] = []
        self.Meta: Dict[str, object] = {}
        self._load()

    def _load(self) -> None:
        if not self.file_path.exists():
            raise FileNotFoundError(f"Cannot find data file: {self.file_path}")

        if self.file_path.suffix.lower() == ".npz":
            self._load_npz()
        else:
            self._load_hdf5()

    def _load_npz(self) -> None:
        archive = np.load(self.file_path)
        signals = archive["signals"]
        wavelengths = archive.get("wavelengths")
        sampling_rate = float(archive.get("sampling_rate", 40e6))
        repetition = int(archive.get("repetition", 1))

        if wavelengths is None:
            wavelengths = np.linspace(680, 900, signals.shape[0] // repetition)

        meta = LoaderMeta(
            wavelengths=np.asarray(wavelengths),
            repetition=repetition,
            detector_count=signals.shape[1],
            sampling_rate=sampling_rate,
        )
        self._frames = _flatten_frames(signals, wavelengths, repetition)
        self.Meta = meta.__dict__

    def _load_hdf5(self) -> None:
        with h5py.File(self.file_path, "r") as handle:
            dataset = _find_first_dataset(handle)
            data = np.array(dataset)

            if data.ndim == 2:
                data = data[np.newaxis, ...]

            wavelengths = _read_wavelengths(handle, data.shape[0])
            repetition = int(handle.attrs.get("RepNum", 1))
            sampling_rate = float(handle.attrs.get("SamplingRate", 40e6))

            meta = LoaderMeta(
                wavelengths=wavelengths,
                repetition=repetition,
                detector_count=data.shape[1],
                sampling_rate=sampling_rate,
            )

        self._frames = _flatten_frames(data, wavelengths, repetition)
        self.Meta = meta.__dict__

    def __len__(self) -> int:  # pragma: no cover - trivial wrapper
        return len(self._frames)

    def __iter__(self) -> Iterator[Frame]:  # pragma: no cover - trivial wrapper
        return iter(self._frames)

    def __call__(self, index: int) -> Frame:
        return self._frames[index]


def _find_first_dataset(handle: h5py.File) -> h5py.Dataset:
    for _, dataset in handle.items():
        if isinstance(dataset, h5py.Dataset):
            return dataset
        if isinstance(dataset, h5py.Group):
            result = _find_first_dataset(dataset)
            if result is not None:
                return result
    raise RuntimeError("No dataset found in HDF5 file")


def _read_wavelengths(handle: h5py.File, fallback_length: int) -> np.ndarray:
    if "wavelength" in handle:
        wavelengths = np.array(handle["wavelength"])
        return wavelengths.reshape(-1)

    for name, dataset in handle.items():
        if name.lower().startswith("wavelength"):
            wavelengths = np.array(dataset)
            return wavelengths.reshape(-1)

    return np.linspace(680, 900, fallback_length)


def _flatten_frames(signals: np.ndarray, wavelengths: np.ndarray, repetition: int) -> List[Frame]:
    frames: List[Frame] = []
    total_frames = signals.shape[0]
    expected = len(wavelengths) * repetition
    if expected == total_frames:
        wavelengths = np.repeat(wavelengths, repetition)

    for idx, signal in enumerate(signals):
        frame_meta = {"Wavelength": float(wavelengths[idx % len(wavelengths)])}
        frames.append(Frame(signal.astype(np.float32), frame_meta))
    return frames


def create_synthetic_dataset(
    wavelengths: np.ndarray,
    detector_count: int = 256,
    samples: int = 512,
    repetition: int = 1,
) -> np.ndarray:
    """Generate a synthetic dataset useful for testing the pipeline."""

    total_frames = len(wavelengths) * repetition
    time = np.linspace(0, 1e-5, samples)
    signals = np.zeros((total_frames, detector_count, samples), dtype=np.float32)

    for frame_idx, wavelength in enumerate(np.repeat(wavelengths, repetition)):
        frequency = 200e3 + 2e3 * np.cos(np.deg2rad(wavelength - wavelengths.mean()))
        envelope = np.exp(-((time - time.mean()) ** 2) / (2 * (1e-6) ** 2))
        phase = 2 * np.pi * frequency * time
        signal = envelope * np.sin(phase)
        signals[frame_idx] = signal
    return signals
