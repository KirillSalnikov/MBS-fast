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
    def test_result_mount_failure_is_not_hidden_by_writable_queue(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);lane=root/'lane';output=root/'other_mount'
            lane.mkdir();output.mkdir()
            with patch.object(guard,'check_storage',side_effect=[
                    dict(ready=True),dict(ready=False,reason='metadata full',errno=28)]) as probe:
                health=guard.check_job_storage(lane,output)
            self.assertFalse(health['ready'])
            self.assertEqual(health['role'],'result')
            self.assertEqual(health['path'],str(output))
            self.assertEqual(probe.call_count,2)
    def test_new_relative_result_checks_existing_parent(self):
        with tempfile.TemporaryDirectory() as d:
            lane=Path(d);(lane/'results').mkdir()
            with patch.object(guard,'check_storage',return_value=dict(ready=True)) as probe:
                self.assertTrue(guard.check_job_storage(lane,'results/new_task')['ready'])
            self.assertEqual(probe.call_args_list[-1].args[0],lane/'results')
    def test_permission_error_is_a_storage_block(self):
        with tempfile.TemporaryDirectory() as d,patch.object(guard,'check_storage',side_effect=PermissionError(13,'denied')):
            health=guard.check_job_storage(d)
            self.assertFalse(health['ready']);self.assertEqual(health['errno'],13)
if __name__=='__main__':unittest.main()
