from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass, field, asdict
from typing import List, Optional

import numpy as np

from .data_models import ReconstructionResult
from .filtering import MSOTPreFilter, PreFilterSettings
from .loader import MSOTSignalLoader, create_synthetic_dataset
from .reconstruction import ReconSettings, ReconSystem, tune_speed_of_sound
from .unmixing import UnmixSettings, UnmixSystem


@dataclass
class PipelineSettings:
    input_path: pathlib.Path
    output_dir: pathlib.Path
    recon: ReconSettings = field(default_factory=ReconSettings)
    filter: PreFilterSettings = field(default_factory=PreFilterSettings)
    unmix: Optional[UnmixSettings] = None
    tune_speed: bool = False

    @classmethod
    def from_json(cls, path: str | pathlib.Path) -> "PipelineSettings":
        raw = json.loads(pathlib.Path(path).read_text())
        recon = ReconSettings(**raw.get("recon", {}))
        filt = PreFilterSettings(**raw.get("filter", {}))
        unmix_raw = raw.get("unmix")
        unmix = UnmixSettings(**unmix_raw) if unmix_raw else None
        return cls(
            input_path=pathlib.Path(raw["input_path"]),
            output_dir=pathlib.Path(raw["output_dir"]),
            recon=recon,
            filter=filt,
            unmix=unmix,
            tune_speed=bool(raw.get("tune_speed", False)),
        )

    def to_json(self, path: str | pathlib.Path) -> None:
        data = asdict(self)
        data["input_path"] = str(self.input_path)
        data["output_dir"] = str(self.output_dir)
        pathlib.Path(path).write_text(json.dumps(data, indent=2))


def run_pipeline(settings: PipelineSettings) -> dict:
    """Execute the preprocessing, reconstruction, and unmixing steps."""

    loader = MSOTSignalLoader(settings.input_path)
    filter_stage = MSOTPreFilter(settings.filter)

    if settings.tune_speed:
        tuned_speed = tune_speed_of_sound(loader, settings.recon.field_of_view, settings.recon.image_size[0])
        settings.recon.speed_of_sound = tuned_speed

    reconstructor = ReconSystem(settings.recon)

    reconstructions: List[ReconstructionResult] = []
    for frame in loader:
        filtered = filter_stage(frame)
        recon_result = reconstructor(filtered)
        reconstructions.append(recon_result)

    wavelengths = np.array([frame.meta.get("Wavelength") for frame in loader])
    output = {
        "reconstructions": np.stack([r.image for r in reconstructions]),
        "wavelengths": wavelengths,
    }

    if settings.unmix:
        unmixer = UnmixSystem(settings.unmix, wavelengths)
        unmix_result = unmixer(reconstructions)
        output["unmixed"] = unmix_result.unmixed_image

    settings.output_dir.mkdir(parents=True, exist_ok=True)
    np.save(settings.output_dir / "reconstructions.npy", output["reconstructions"])
    np.save(settings.output_dir / "wavelengths.npy", wavelengths)
    if "unmixed" in output:
        np.save(settings.output_dir / "unmixed.npy", output["unmixed"])

    with open(settings.output_dir / "pipeline_settings.json", "w") as handle:
        handle.write(json.dumps(asdict(settings), default=str, indent=2))

    return output


def demo(output_dir: str | pathlib.Path, repetition: int = 1) -> dict:
    """Run the pipeline using synthetic data to demonstrate usage."""

    wavelengths = np.array([700, 730, 760, 800, 850], dtype=np.float32)
    signals = create_synthetic_dataset(wavelengths, repetition=repetition)
    input_path = pathlib.Path(output_dir) / "synthetic.npz"
    np.savez(input_path, signals=signals, wavelengths=wavelengths, repetition=repetition)

    settings = PipelineSettings(
        input_path=input_path,
        output_dir=pathlib.Path(output_dir),
        recon=ReconSettings(),
        filter=PreFilterSettings(),
        unmix=UnmixSettings(
            endmember_names=["Hb", "HbO2"],
            spectra={
                "Hb": np.linspace(0.6, 0.9, len(wavelengths)),
                "HbO2": np.linspace(1.0, 0.7, len(wavelengths)),
            },
        ),
        tune_speed=False,
    )
    return run_pipeline(settings)
