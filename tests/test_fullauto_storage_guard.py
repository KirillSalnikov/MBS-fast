import importlib.util
import os
import subprocess
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

class BtrfsAdmissionTests(unittest.TestCase):
    fsid = 'fbd468ae-3397-41ec-b376-11cdd1b6db2a'
    # Read-only GPU1 sysfs snapshot, 2026-11-05. Global reserve is already in
    # bytes_may_use; df simultaneously reported more than63GB available.
    gpu1 = dict(total_bytes=48398073856, bytes_used=47856336896,
                bytes_may_use=538312704, bytes_reserved=0, bytes_pinned=0,
                bytes_readonly=0, bytes_zone_unusable=0, disk_total=48398073856,
                flags=4)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.allocation = self.root / self.fsid / 'allocation'
        (self.allocation / 'metadata').mkdir(parents=True)
        for key, value in self.gpu1.items():
            (self.allocation / 'metadata' / key).write_text(str(value))
        for name in ('global_rsv_size', 'global_rsv_reserved'):
            (self.allocation / name).write_text('536870912')
        self.sysfs = patch.object(guard, 'BTRFS_SYSFS', self.root)
        self.identity = patch.object(guard, '_btrfs_identity', return_value=(self.fsid, 1))
        self.sysfs.start(); self.identity.start()
        self.addCleanup(self.sysfs.stop); self.addCleanup(self.identity.stop)

    def change(self, name, value):
        (self.allocation / 'metadata' / name).write_text(str(value))

    def probe(self, unallocated=1048576):
        with patch.object(guard, '_btrfs_usage', return_value=dict(device_unallocated_bytes=unallocated)):
            return guard.btrfs_headroom('/')

    def test_actual_gpu1_blocks_before_any_temp_write_despite_df_space(self):
        st = SimpleNamespace(f_bavail=63058903040, f_frsize=1, f_files=0, f_favail=0)
        with patch.object(guard.os, 'statvfs', return_value=st), \
             patch.object(guard, 'filesystem_type', return_value='btrfs'), \
             patch.object(guard, '_btrfs_usage', return_value=dict(device_unallocated_bytes=1048576)), \
             patch.object(guard.tempfile, 'mkstemp') as create:
            health = guard.check_storage('/')
        self.assertFalse(health['ready']); create.assert_not_called()
        self.assertEqual(health['btrfs']['metadata_headroom_bytes'], 3424256)
        self.assertEqual(health['btrfs']['global_rsv_reserved'], 536870912)
        self.assertIn('metadata reserve', health['reason'])

    def test_reserved_readonly_pinned_and_zoned_bytes_all_count(self):
        for name, value in [('bytes_reserved', 100), ('bytes_pinned', 200),
                            ('bytes_readonly', 300), ('bytes_zone_unusable', 400)]:
            self.change(name, value)
        self.assertEqual(self.probe()['metadata_headroom_bytes'], 3424256-1000)

    def test_negative_headroom_is_reported_not_clamped(self):
        self.change('bytes_reserved', 8*1024**2)
        health = self.probe()
        self.assertFalse(health['ready'])
        self.assertEqual(health['metadata_headroom_bytes'], 3424256-8*1024**2)

    def test_one_free_chunk_cannot_cover_overcommitted_metadata(self):
        self.change('bytes_reserved', 3*1024**3)
        with patch.object(guard, '_btrfs_usage', return_value=dict(device_unallocated_bytes=1024**3)) as usage:
            health = guard.btrfs_headroom('/')
        self.assertFalse(health['ready'])
        self.assertIn('insufficient metadata', health['reason'])
        self.assertEqual(health['allocation_margin_bytes'],3*1024**3-3424256+1024**3)
        usage.assert_called_once()

    def test_negative_slack_can_use_known_single_device_backing(self):
        self.change('bytes_reserved',3*1024**3)
        deficit=3*1024**3-3424256;required=deficit+1024**3
        self.assertFalse(self.probe(required-1)['ready'])
        health=self.probe(required)
        self.assertTrue(health['ready'])
        self.assertEqual(health['metadata_headroom_bytes'],-deficit)
        self.assertEqual(health['metadata_reservation_deficit_bytes'],deficit)
        self.assertEqual(health['physical_metadata_ratio'],1)
        self.assertIn('not guaranteed',health['warning'])

    def test_negative_dup_slack_needs_two_physical_copies(self):
        self.change('bytes_reserved',512*1024**2)
        self.change('disk_total',self.gpu1['total_bytes']*2)
        logical=512*1024**2-3424256+1024**3
        self.assertFalse(self.probe(logical)['ready'])
        self.assertFalse(self.probe(2*logical-1)['ready'])
        self.assertTrue(self.probe(2*logical)['ready'])

    def test_negative_slack_unknown_or_multidevice_placement_blocks(self):
        self.change('bytes_reserved',512*1024**2)
        with patch.object(guard,'_btrfs_identity',return_value=(self.fsid,2)):
            self.assertFalse(self.probe(100*1024**3)['ready'])
        self.change('disk_total',self.gpu1['total_bytes']*3)
        self.assertFalse(self.probe(100*1024**3)['ready'])
        self.change('disk_total',self.gpu1['total_bytes']+1)
        self.assertFalse(self.probe(100*1024**3)['ready'])
        self.change('disk_total',self.gpu1['total_bytes'])
        self.change('flags',12)
        self.assertFalse(self.probe(100*1024**3)['ready'])

    def test_gpu2_recorded_overcommit_with_unallocated_capacity(self):
        values=dict(total_bytes=82686509056,bytes_used=80124772352,
                    bytes_may_use=11941183488,bytes_reserved=345997312,
                    bytes_pinned=87851008,bytes_readonly=65536,
                    bytes_zone_unusable=0,disk_total=82686509056,flags=4)
        for name,value in values.items():self.change(name,value)
        health=self.probe(126253793280)
        self.assertTrue(health['ready'])
        self.assertEqual(health['metadata_headroom_bytes'],-9813360640)
        self.assertEqual(health['allocation_margin_bytes'],9813360640+1024**3)

    def test_negative_slack_missing_usage_fails_closed(self):
        self.change('bytes_reserved',512*1024**2)
        with patch.object(guard,'_btrfs_usage',side_effect=RuntimeError('unknown capacity')):
            health=guard.btrfs_headroom('/')
        self.assertFalse(health['ready']);self.assertIn('information unavailable',health['reason'])

    def test_healthy_existing_metadata_does_not_spawn_usage(self):
        self.change('bytes_used', self.gpu1['bytes_used']-256*1024**2)
        with patch.object(guard, '_btrfs_usage') as usage:
            self.assertTrue(guard.btrfs_headroom('/')['ready'])
        usage.assert_not_called()

    def test_single_device_allocator_margin_is_explicit_not_guaranteed(self):
        health = self.probe(1024**3)
        self.assertTrue(health['ready'])
        self.assertIn('not guaranteed', health['warning'])
        self.assertEqual(health['allocation_margin_bytes'], 1024**3)

    def test_dup_profile_needs_physical_room_for_both_copies(self):
        self.change('disk_total', self.gpu1['total_bytes']*2)
        self.assertFalse(self.probe(1024**3)['ready'])
        self.assertTrue(self.probe(2*1024**3)['ready'])

    def test_multidevice_low_headroom_does_not_infer_placement_from_total(self):
        with patch.object(guard, '_btrfs_identity', return_value=(self.fsid, 2)):
            health = self.probe(100*1024**3)
        self.assertFalse(health['ready']); self.assertIn('placement unknown', health['reason'])

    def test_mixed_group_low_headroom_is_not_misread_as_metadata_only(self):
        self.change('flags', 5)
        self.assertIn('placement unknown', self.probe(100*1024**3)['reason'])

    def test_missing_reservations_fail_closed_for_known_btrfs(self):
        (self.allocation / 'metadata/bytes_may_use').unlink()
        health = self.probe()
        self.assertFalse(health['ready']); self.assertIn('information unavailable', health['reason'])

    def test_missing_optional_zoned_counter_on_older_kernel_is_zero(self):
        (self.allocation / 'metadata/bytes_zone_unusable').unlink()
        self.assertEqual(self.probe()['metadata_headroom_bytes'], 3424256)

    def test_usage_timeout_is_explicit_unknown_and_blocks(self):
        with patch.object(guard, '_btrfs_usage', side_effect=subprocess.TimeoutExpired(['btrfs'], 3)):
            health = guard.btrfs_headroom('/')
        self.assertFalse(health['ready']); self.assertIn('information unavailable', health['reason'])

    def test_healthy_ext4_does_not_query_btrfs(self):
        with tempfile.TemporaryDirectory() as d, \
             patch.object(guard.os, 'statvfs', return_value=SimpleNamespace(f_bavail=10*1024**3,f_frsize=1,f_files=1000000,f_favail=900000)), \
             patch.object(guard, 'filesystem_type', return_value='ext4'), \
             patch.object(guard, 'btrfs_headroom') as btrfs:
            self.assertTrue(guard.check_storage(d)['ready'])
            self.assertEqual(list(Path(d).iterdir()), [])
        btrfs.assert_not_called()

    def test_usage_parser_uses_C_locale_raw_bytes_and_timeout(self):
        output = 'Overall:\n    Device unallocated:\t 1048576\nMetadata,single: Size:48398073856, Used:47856336896 (98.88%)\n'
        with patch.object(guard.os.path, 'isfile', return_value=True), \
             patch.object(guard.os, 'access', return_value=True), \
             patch.object(guard.subprocess, 'run', return_value=SimpleNamespace(returncode=0,stdout=output,stderr='WARNING: run as root for detailed info')) as command:
            health = guard._btrfs_usage('/path with spaces')
        self.assertEqual(health['device_unallocated_bytes'], 1048576)
        self.assertEqual(command.call_args[1]['timeout'], 3)
        self.assertEqual(command.call_args[1]['env']['LC_ALL'], 'C')
        self.assertEqual(command.call_args[0][0][-1], '/path with spaces')
        self.assertNotIn('sudo', command.call_args[0][0])

    def test_fs_info_binds_sysfs_to_actual_filesystem_identity(self):
        self.identity.stop()
        def ioctl(fd, request, result, mutate):
            self.assertEqual(request, 0x8400941f)
            result[8:16] = (1).to_bytes(8, 'little')
            result[16:32] = __import__('uuid').UUID(self.fsid).bytes
        with tempfile.TemporaryDirectory() as d, patch.object(guard.fcntl, 'ioctl', side_effect=ioctl):
            self.assertEqual(guard._btrfs_identity(d), (self.fsid, 1))

    def test_mountinfo_respects_deep_bind_and_escaped_path(self):
        text = '1 0 0:1 / / rw - btrfs /dev/sdb2 rw\n2 1 8:1 / /home/data\\040copy rw - ext4 /dev/sda1 rw\n'
        with patch.object(guard.Path, 'read_text', return_value=text):
            self.assertEqual(guard.filesystem_type('/home/data copy/subdir'), 'ext4')
            self.assertEqual(guard.filesystem_type('/home/data copycat'), 'btrfs')


class PressureConfirmationTests(unittest.TestCase):
    def pressure(self, slack=-458752):
        return dict(ready=False, reason='Btrfs metadata reservations exceed allocated space',
                    role='result', path='/result', btrfs=dict(
                        ready=False, failure_kind='metadata_pressure', metadata_headroom_bytes=slack))

    def test_transient_low_then_ready_keeps_native_admitted(self):
        with patch.object(guard, 'check_job_storage', side_effect=[self.pressure(), dict(ready=True, reason='verified')]) as probe, \
             patch.object(guard.time, 'sleep') as sleep:
            health = guard.check_job_storage_confirmed('/lane', '/result')
        self.assertTrue(health['ready']); self.assertEqual(probe.call_count, 2)
        sleep.assert_called_once_with(2)
        self.assertEqual(health['storage_confirmation']['samples'][0]['btrfs']['metadata_headroom_bytes'], -458752)
        self.assertFalse(health['storage_confirmation']['persistent_allocator_pressure'])

    def test_persistent_pressure_is_bounded_and_fails_closed(self):
        with patch.object(guard, 'check_job_storage', return_value=self.pressure()) as probe, \
             patch.object(guard.time, 'sleep') as sleep:
            health = guard.check_job_storage_confirmed('/lane', '/result')
        self.assertFalse(health['ready']); self.assertEqual(probe.call_count, 4)
        self.assertEqual(sleep.call_count, 3)
        self.assertEqual(health['storage_confirmation']['wait_seconds'], 6)
        self.assertTrue(health['storage_confirmation']['persistent_allocator_pressure'])

    def test_actual_io_failure_is_immediate_even_after_pressure(self):
        error = dict(ready=False, errno=28, reason='No space left on device')
        for sequence, expected_calls in [([error], 1), ([self.pressure(), error], 2)]:
            with patch.object(guard, 'check_job_storage', side_effect=sequence) as probe, \
                 patch.object(guard.time, 'sleep') as sleep:
                health = guard.check_job_storage_confirmed('/lane', '/result')
            self.assertEqual(health['errno'], 28); self.assertEqual(probe.call_count, expected_calls)
            self.assertEqual(sleep.call_count, expected_calls-1)

    def test_missing_allocator_info_or_byte_inode_reserve_does_not_wait(self):
        for reason in ['inode reserve', 'byte reserve', 'Btrfs allocation information unavailable']:
            with patch.object(guard, 'check_job_storage', return_value=dict(ready=False, reason=reason)) as probe, \
                 patch.object(guard.time, 'sleep') as sleep:
                self.assertFalse(guard.check_job_storage_confirmed('/lane')['ready'])
            self.assertEqual(probe.call_count, 1); sleep.assert_not_called()

    def pinned(self):
        health = self.pressure(57622528)
        health['btrfs']['metadata'] = dict(bytes_pinned=175276032)
        return health

    def test_recorded_pinned_pressure_sync_rechecks_all_original_checks(self):
        with patch.object(guard, 'check_job_storage', side_effect=[self.pinned()]*4+[dict(ready=True)]) as probe, \
             patch.object(guard.time, 'sleep'), \
             patch.object(guard, 'drain_btrfs_transaction', return_value=dict(success=True)) as sync:
            health = guard.check_job_storage_confirmed('/lane','/result')
        self.assertTrue(health['ready']); self.assertEqual(probe.call_count,5)
        sync.assert_called_once_with('/result')
        self.assertEqual(health['storage_confirmation']['attempts'],5)
        self.assertFalse(health['storage_confirmation']['persistent_allocator_pressure'])

    def test_successful_sync_is_not_admission_if_pressure_remains(self):
        with patch.object(guard, 'check_job_storage', return_value=self.pinned()) as probe, \
             patch.object(guard.time, 'sleep'), \
             patch.object(guard, 'drain_btrfs_transaction', return_value=dict(success=True)) as sync:
            health = guard.check_job_storage_confirmed('/lane','/result')
        self.assertFalse(health['ready']); self.assertEqual(probe.call_count,5)
        sync.assert_called_once()

    def test_sync_failure_does_not_retry_or_admit(self):
        with patch.object(guard, 'check_job_storage', return_value=self.pinned()) as probe, \
             patch.object(guard.time, 'sleep'), \
             patch.object(guard, 'drain_btrfs_transaction', return_value=dict(success=False,error='I/O error')) as sync:
            health=guard.check_job_storage_confirmed('/lane','/result')
        self.assertFalse(health['ready']);self.assertEqual(probe.call_count,4)
        self.assertIn('I/O error',health['reason']);sync.assert_called_once()

    def test_real_io_error_and_unpinned_pressure_never_sync(self):
        for health in (dict(self.pinned(),errno=5),self.pressure()):
            with patch.object(guard,'check_job_storage',return_value=health), \
                 patch.object(guard.time,'sleep'),patch.object(guard,'drain_btrfs_transaction')as sync:
                self.assertFalse(guard.check_job_storage_confirmed('/lane','/result')['ready'])
            sync.assert_not_called()

    def test_post_sync_io_error_is_preserved(self):
        with patch.object(guard,'check_job_storage',side_effect=[self.pinned()]*4+[dict(ready=False,errno=5,reason='EIO')]), \
             patch.object(guard.time,'sleep'),patch.object(guard,'drain_btrfs_transaction',return_value=dict(success=True)):
            health=guard.check_job_storage_confirmed('/lane','/result')
        self.assertFalse(health['ready']);self.assertEqual(health['errno'],5)

    def test_sync_stderr_receipt_reports_actual_post_admission(self):
        after=dict(ready=True,checks=[dict(path='/result',btrfs=dict(metadata_headroom_bytes=128*1024**2))])
        with patch.object(guard,'check_job_storage',side_effect=[self.pinned()]*4+[after]), \
             patch.object(guard.time,'sleep'),patch.object(guard,'drain_btrfs_transaction',return_value=dict(success=True,seconds=.25)), \
             patch.object(guard.os,'write')as output:
            health=guard.check_job_storage_confirmed('/lane','/result')
        self.assertTrue(health['ready']);output.assert_called_once()
        fd,line=output.call_args[0];self.assertEqual(fd,2)
        payload=__import__('json').loads(line.decode()[len('STORAGE transaction_sync '):])
        self.assertEqual(payload['path'],'/result');self.assertTrue(payload['success']);self.assertTrue(payload['post_ready'])
        self.assertEqual(payload['post_metadata_headroom_bytes'],128*1024**2)

    def test_stderr_failure_preserves_failed_sync_hold(self):
        with patch.object(guard,'check_job_storage',return_value=self.pinned()), \
             patch.object(guard.time,'sleep'),patch.object(guard,'drain_btrfs_transaction',return_value=dict(success=False,error='timeout')), \
             patch.object(guard.os,'write',side_effect=OSError(28,'log full')):
            health=guard.check_job_storage_confirmed('/lane','/result')
        self.assertFalse(health['ready']);self.assertIn('timeout',health['reason'])

class TransactionSyncTests(unittest.TestCase):
    def test_success_targets_one_filesystem_without_shell_or_sudo(self):
        with patch.object(guard.os.path,'isfile',return_value=True),patch.object(guard.os,'access',return_value=True), \
             patch.object(guard.subprocess,'Popen')as spawn:
            spawn.return_value.wait.return_value=0
            health=guard.drain_btrfs_transaction('/result with spaces')
        self.assertTrue(health['success'])
        self.assertEqual(spawn.call_args[0][0],['/usr/bin/sync','-f','/result with spaces'])
        spawn.return_value.wait.assert_called_once_with(timeout=10)

    def test_hung_child_has_bounded_reap_and_reports_failure(self):
        with patch.object(guard.os.path,'isfile',return_value=True),patch.object(guard.os,'access',return_value=True), \
             patch.object(guard.subprocess,'Popen')as spawn:
            child=spawn.return_value;child.poll.return_value=None
            child.wait.side_effect=subprocess.TimeoutExpired(['sync'],10)
            health=guard.drain_btrfs_transaction('/result')
        self.assertFalse(health['success']);self.assertTrue(health['child_pending'])
        child.kill.assert_called_once()
        self.assertEqual([c.kwargs['timeout']for c in child.wait.call_args_list],[10,1])

    def test_nonzero_sync_is_failure(self):
        with patch.object(guard.os.path,'isfile',return_value=True),patch.object(guard.os,'access',return_value=True), \
             patch.object(guard.subprocess,'Popen')as spawn:
            spawn.return_value.wait.return_value=1;spawn.return_value.poll.return_value=1
            health=guard.drain_btrfs_transaction('/result')
        self.assertFalse(health['success']);self.assertIn('exit 1',health['error'])

if __name__=='__main__':unittest.main()
