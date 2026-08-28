"""Reproduce the manuscript simulations, validation analyses, and figures."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent
PYTHON_CODE = ROOT / "Python_code"
FIGURE_DIR = ROOT / "outputs" / "manuscript_figures"
sys.path.insert(0, str(PYTHON_CODE))


def run_python_analyses() -> None:
    """Run the conceptual figure, onset, recovery, and sparse-trajectory analyses."""
    import export_figure1_concept
    import sim_psychopathology
    import recovery_bistable
    import clinical_trajectory_comparison

    export_figure1_concept.main()
    sim_psychopathology.run_all()
    recovery_bistable.run_all()
    clinical_trajectory_comparison.run_all()


def run_r_analyses() -> None:
    """Run the intensive-data validation and regenerate the Wichers figure."""
    subprocess.run(
        ["Rscript", str(ROOT / "VALIDATION" / "run_validation.R")],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        ["Rscript", str(ROOT / "VALIDATION" / "generate_wichers_figure.R")],
        cwd=ROOT,
        check=True,
    )


def verify_outputs() -> None:
    """Require all six canonical manuscript figures in raster and vector form."""
    expected = [
        FIGURE_DIR / f"Fig_{number}{suffix}"
        for number in range(1, 6)
        for suffix in (".png", ".pdf")
    ] + [
        FIGURE_DIR / f"Fig_6_Wichers{suffix}"
        for suffix in (".png", ".pdf")
    ]
    missing = [path.relative_to(ROOT) for path in expected if not path.is_file()]
    if missing:
        raise RuntimeError("Missing expected outputs: " + ", ".join(map(str, missing)))


def main() -> int:
    """Run the complete public reproducibility workflow."""
    started = time.time()
    run_python_analyses()
    run_r_analyses()
    verify_outputs()
    elapsed = (time.time() - started) / 60.0
    print(f"Complete reproducibility workflow finished in {elapsed:.1f} minutes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
