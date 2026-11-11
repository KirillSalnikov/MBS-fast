import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('archive', Path(__file__).resolve().parents[1] / 'scripts/archive_fullauto_results.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class ArchiveTest(unittest.TestCase):
    def fixture(self, root):
        (root / '.mbs-result-layout-v1').write_text('mbs_compact_result_v1\n')
        (root / 'work').mkdir()
        (root / 'work/fullauto_checkpoint.bin').write_bytes(b'recovery data')
        (root / 'results').mkdir()
        for i in range(5):
            p = root / 'results' / ('snapshot_' + str(i))
            p.mkdir()
            (p / 'manifest.fnv').write_text('fixture\n')
            (p / 'mueller_bins.dat').write_text('matrix ' + str(i))
        (root / 'latest').symlink_to('results/snapshot_4')

    def test_roundtrip(self):
        for compression in ('gzip', 'zstd'):
            if compression == 'zstd' and not a.shutil.which('zstd'):
                continue
            with tempfile.TemporaryDirectory() as t:
                root = Path(t)
                self.fixture(root)
                dry = a.archive_results(root)
                self.assertEqual(len(dry['selected']), 3)
                r = a.archive_results(root, apply=True, compression=compression)
                self.assertEqual(len(list((root / 'results').iterdir())), 2)
                self.assertEqual((root / 'latest/mueller_bins.dat').read_text(), 'matrix 4')
                self.assertEqual((root / 'work/fullauto_checkpoint.bin').read_bytes(), b'recovery data')
                receipt = json.loads(Path(r['archive'] + '.json').read_text())
                with a.reader(Path(r['archive']), compression) as tar:
                    found = {m.name for m in tar if m.isfile()}
                self.assertEqual(found, set(receipt['files']))
                self.assertFalse(a.archive_results(root, apply=True, compression=compression)['applied'])

    def test_changed_source_retained(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.fixture(root)
            real = a.inventory
            count = [0]
            def changed(*args):
                count[0] += 1
                result = real(*args)
                return result if count[0] == 1 else {}
            with patch.object(a, 'inventory', changed):
                with self.assertRaises(ValueError):
                    a.archive_results(root, apply=True, compression='gzip')
            self.assertEqual(len(list((root / 'results').iterdir())), 5)

    def test_symlink_refused(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.fixture(root)
            (root / 'results/snapshot_0/unexpected').symlink_to(root / 'work/fullauto_checkpoint.bin')
            a.os.utime(str(root / 'results/snapshot_0'), (1, 1))
            with self.assertRaises(ValueError):
                a.archive_results(root, apply=True, compression='gzip')
            self.assertTrue((root / 'work/fullauto_checkpoint.bin').exists())

    def test_active_lock_refused(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.fixture(root)
            with (root / 'work/controller.lock').open('a') as f:
                a.fcntl.flock(f, a.fcntl.LOCK_EX | a.fcntl.LOCK_NB)
                with self.assertRaises(BlockingIOError):
                    a.archive_results(root, apply=True, compression='gzip')


if __name__ == '__main__':
    unittest.main()
