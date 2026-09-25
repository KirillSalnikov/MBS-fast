#!/usr/bin/env python3
"""CPU integral-only regression; uses only the Python standard library."""
import argparse
import csv
import math
import pathlib
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=pathlib.Path, default=ROOT / "bin/mbs_po")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="mbs-integrals-test-") as tmp:
        work = pathlib.Path(tmp)
        shape = work / "cube.particle"
        # Compare identical FULL-domain orientation sets, not cube symmetry sets.
        shape.write_text((ROOT / "examples/cube.particle").read_text().replace("90 90", "180 360"))
        common = [str(args.binary.resolve()), "--method", "po", "--backend", "cpu",
                  "--geometry", "convex", "--particle-file", str(shape),
                  "--wavelength-um", "1.064", "--hammersley", "64",
                  "--max-reflections", "8", "--beam-cutoff-jones", "0.001",
                  "--beam-cutoff-area", "0.002", "--trace-cutoff-importance", "0.0001",
                  "--trace-max-beams", "20000", "--trace-limit-retries", "0", "--close"]

        def run(name, options, error=False):
            output = work / name
            result = subprocess.run(common + options + ["--output", str(output)],
                                    capture_output=True, text=True, timeout=120)
            if error:
                assert result.returncode != 0, name
                assert "--integrals-only" in result.stdout + result.stderr, result.stdout + result.stderr
                return
            assert result.returncode == 0, result.stdout + result.stderr
            with next(output.glob("*integrals.tsv")).open() as stream:
                return list(csv.DictReader(stream, delimiter="\t"))

        for index in ["0.000038", "0.0723"]:
            physics = ["--refractive-index", "1.5563", index, "--resize-dmax-um", "6"]
            fast = run("fast" + index, physics + ["--integrals-only", "--threads", "4"])[0]
            serial = run("serial" + index, physics + ["--integrals-only", "--threads", "1"])[0]
            assert {k: v for k, v in fast.items() if k != "seconds"} == {
                k: v for k, v in serial.items() if k != "seconds"}
            full = run("full" + index, physics + ["--threads", "4", "--scattering-grid", "0", "180", "12", "90"])[0]
            assert math.isclose(float(fast["Cext_OT"]), float(full["Cext"]), rel_tol=1e-12)
            print("PASS: identical forward Cext and thread-independent reduction, ni=" + index)

        physics = ["--refractive-index", "1.5563", "0.000038", "--integrals-only", "--threads", "4"]
        shared = run("shared", physics + ["--dmax-grid", "6", "12", "2"])
        single = run("single", physics + ["--resize-dmax-um", "12"])[0]
        assert len(shared) == 2
        assert math.isclose(float(shared[-1]["Dmax_um"]), 12, rel_tol=1e-6)
        assert math.isclose(float(shared[-1]["Cext_OT"]), float(single["Cext_OT"]), rel_tol=1e-4)
        print("PASS: shared size geometry versus independent size")
        for i, extra in enumerate([["--symmetry", "2", "2"], ["--mirror-gamma"],
                                   ["--incoherent"], ["--no-shadow-beam"],
                                   ["--ot-phase-average"]]):
            run("reject" + str(i), physics + extra, error=True)
        print("PASS: incompatible integral-only options rejected")


if __name__ == "__main__":
    main()
