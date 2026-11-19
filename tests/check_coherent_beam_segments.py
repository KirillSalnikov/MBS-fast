#!/usr/bin/env python3
"""Compare dense coherent CUDA batches to CPU FP64 on a joined bullet rosette."""
import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile


def calculate(binary, backend, batch, root):
    output = root / (backend + str(batch))
    command = [
        str(binary), '--method', 'po', '--backend', backend,
        '--particle', '3', '1.887064518871375', '3.4463361198126723',
        '--rosette-joined', '--symmetry', '2', '2', '--wavelength-um', '.532',
        '--ri', '1.3116', '1.4898e-9', '--haar-alpha', '--sobol-seed', '16', '2003',
        '--scattering-grid', '0', '180', '16', '8', '--max-reflections', '8',
        '--threads', '1', '--no-gpu-trace-prefilter', '--po-trace-backend', 'cpu',
        '--orientation-chunk', str(batch), '--close', '--output', str(output),
    ]
    log = subprocess.check_output(command, cwd=str(root), stderr=subprocess.STDOUT,
                                  universal_newlines=True)
    if backend == 'cuda':
        assert 'Phase 2 (diffraction, CUDA)' in log, 'CUDA diffraction was not used'
    beams = int(re.search(r'Candidate output beams:\s*(\d+)', log).group(1))
    assert beams > 200000
    with (output / (output.name + '.dat')).open() as stream:
        next(stream)
        rows = [list(map(float, line.split())) for line in stream if line.strip()]
    assert len(rows) == 9 and rows[0][0] == 0 and rows[-1][0] == 180
    assert all(len(x) == 18 and x[2] > 0 and all(map(math.isfinite, x)) for x in rows)
    return beams, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    binary = parser.parse_args().binary.resolve()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        expected_beams, reference = calculate(binary, 'cpu', 16, root)
        for batch in (1, 16):
            beams, rows = calculate(binary, 'cuda', batch, root)
            assert beams == expected_beams and all(x[:2] == y[:2] for x, y in zip(reference, rows))
            raw = max(abs(a - b) / x[2] for x, y in zip(reference, rows)
                      for a, b in zip(x[2:], y[2:]))
            ratio = max(abs(a / x[2] - b / y[2]) for x, y in zip(reference, rows)
                        for a, b in zip(x[3:], y[3:]))
            print('batch=%d max errors: raw/M11 %.8g, normalized %.8g' % (batch, raw, ratio))
            assert max(raw, ratio) < 1e-4


if __name__ == '__main__':
    main()
