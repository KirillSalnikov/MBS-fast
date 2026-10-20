#!/usr/bin/env python3
"""Metadata decoder must retain legacy support and reject damaged compact data."""
import tempfile
import unittest
from pathlib import Path
from fullauto_sampling import fnv, sampling_text


class SamplingFormatTest(unittest.TestCase):
    def test_legacy_and_compact(self):
        data = b"signature identity 128 0 11 8 4 64 0 0 0 0 0 0 auto weights execution\n170\n180\n"
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name)
            (folder / "sampling.index").write_bytes(data)
            self.assertEqual(sampling_text(folder), data.decode())
            (folder / "sampling.index.fnv").write_text(fnv(data))
            self.assertEqual(sampling_text(folder), data.decode())
            (folder / "sampling.index.fnv").write_text("bad")
            with self.assertRaises(ValueError):
                sampling_text(folder)
            (folder / "sampling.index").unlink()
            marker = b"identity 1 datahash none\nMBS_SAMPLING_V1 " + fnv(data).encode() + b" " + str(len(data)).encode() + b"\n" + data
            (folder / "complete.txt").write_bytes(marker)
            self.assertEqual(sampling_text(folder), data.decode())
            diagnostics = b"1 0 0 0.1 0.2 0.3 0 0 0 1\n"
            framed = marker + b"MBS_DIAGNOSTICS_V1 " + fnv(diagnostics).encode() + b"\n" + diagnostics
            (folder / "complete.txt").write_bytes(framed)
            self.assertEqual(sampling_text(folder), data.decode())
            for bad in [marker[:-1], marker.replace(b"170\n", b"171\n"),
                        marker.replace(b"MBS_SAMPLING_V1", b"MBS_SAMPLING_V2"),
                        b"identity 1 datahash none\n", marker.replace(b"\n170\n180\n", b"\n170\n180\nextra\n"),
                        framed.replace(b"0.3", b"0.4")]:
                (folder / "complete.txt").write_bytes(bad)
                with self.assertRaises(ValueError):
                    sampling_text(folder)


if __name__ == "__main__":
    unittest.main()
