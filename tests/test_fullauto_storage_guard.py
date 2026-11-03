import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

p=Path(__file__).resolve().parents[1]/'scripts/fullauto_storage.py'
spec=importlib.util.spec_from_file_location('guard',p);guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)

class StorageTests(unittest.TestCase):
    def stat(self,files,free):
        return SimpleNamespace(f_bavail=1024**3,f_frsize=4096,f_files=files,f_favail=free)
    def test_btrfs_zero_total_inodes_is_not_full(self):
        with tempfile.TemporaryDirectory()as d,patch.object(guard.os,'statvfs',return_value=self.stat(0,0)):
            self.assertTrue(guard.check_storage(d)['ready'])
            self.assertEqual(list(Path(d).iterdir()),[])
    def test_ext4_inode_reserve_blocks_before_writes(self):
        with patch.object(guard.os,'statvfs',return_value=self.stat(1000000,4)),patch.object(guard.tempfile,'mkstemp')as create:
            self.assertEqual(guard.check_storage('/unused')['reason'],'inode reserve');create.assert_not_called()
    def test_metadata_exhaustion_despite_free_bytes(self):
        with tempfile.TemporaryDirectory()as d,patch.object(guard.os,'statvfs',return_value=self.stat(0,0)),patch.object(guard.os,'replace',side_effect=OSError(28,'metadata full')):
            result=guard.check_storage(d);self.assertFalse(result['ready']);self.assertEqual(result['errno'],28)
            self.assertEqual(list(Path(d).iterdir()),[])
    def test_runtime_error_survives_failed_lane(self):
        with tempfile.TemporaryDirectory()as d,patch.dict(os.environ,{'XDG_RUNTIME_DIR':d}):
            guard.record_runtime('/unwritable-lane',{'status':'storage_blocked'})
            self.assertTrue(guard.runtime_path('/unwritable-lane').exists())
if __name__=='__main__':unittest.main()
