"""Example command-line workflow for the Python MSOT toolkit."""

import argparse
import pathlib

from msot_toolkit import PipelineSettings, UnmixSettings, run_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=pathlib.Path, help="Path to .npz or .msot file")
    parser.add_argument("output", type=pathlib.Path, help="Directory for pipeline outputs")
    parser.add_argument(
        "--tune-speed", action="store_true", help="Auto-tune the speed of sound before reconstruction"
    )
    args = parser.parse_args()

    unmix_settings = UnmixSettings(
        endmember_names=["Hb", "HbO2"],
        spectra={
            "Hb": [0.9, 0.8, 0.7, 0.6],
            "HbO2": [0.6, 0.7, 0.8, 0.9],
        },
    )

    settings = PipelineSettings(
        input_path=args.input,
        output_dir=args.output,
        unmix=unmix_settings,
        tune_speed=args.tune_speed,
    )

    run_pipeline(settings)


if __name__ == "__main__":
    main()
