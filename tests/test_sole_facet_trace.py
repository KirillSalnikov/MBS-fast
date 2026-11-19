"""Compare deferred facet ordering with the fully sorted diagnostic trace."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def run(binary, folder, options, diagnostics):
    command = [str(binary), "--method", "po", "--backend", "cpu",
               "--particle", "3", "2.9585957810762125", "4.57505088811751",
               "--rosette-joined", "--symmetry", "2", "2", "--wavelength-um",
               "0.532", "--ri", "1.3116", "1.4898e-9", "--haar-alpha",
               "--scattering-grid", "0", "180", "4", "8", "--threads", "1",
               "--output", str(folder), "--close", *options]
    if diagnostics:
        command.append("--trace-prefilter-stats")
    result = subprocess.run(command, cwd=folder.parent, check=True,
                            env=dict(os.environ, OMP_THREAD_LIMIT="1"),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    beams = re.search(r"Candidate output beams:\s*(\d+)", result.stdout)
    assert beams and int(beams.group(1)) > 0
    return (beams.group(1),
            (folder / (folder.name + ".dat")).read_bytes(),
            (folder / (folder.name + "_integrals.tsv")).read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    binary = parser.parse_args().binary.resolve()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for depth, seed in [(8, 2003), (22, 1700006195)]:
            options = ["--max-reflections", str(depth), "--sobol-seed", "4", str(seed)]
            reference = run(binary, root / ("sorted" + str(depth)), options, True)
            result = run(binary, root / ("deferred" + str(depth)), options, False)
            assert reference == result, "beam count, Mueller matrix or integrals changed"
            print("n={} passed with {} beams".format(depth, result[0]))


if __name__ == "__main__":
    main()
