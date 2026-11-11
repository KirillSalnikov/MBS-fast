"""Bounded allocator confirmation and real owned-process hold/resume regression."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0,str(SCRIPTS))
import run_fullauto_campaign as campaign


class RuntimeStorageTests(unittest.TestCase):
    def test_hold_is_scoped_and_requires_explicit_transition(self):
        record = dict(status='storage_blocked', current_job='first', updated_unix=100,
                      storage=dict(reason='allocator pressure'))
        self.assertIsNotNone(campaign.storage_hold_reason(dict(status='storage_blocked'), 'first', None))
        self.assertIsNotNone(campaign.storage_hold_reason(dict(status='running'), 'first', record))
        self.assertIsNotNone(campaign.storage_hold_reason(dict(status='pending_resume'), 'first', record))
        self.assertIsNotNone(campaign.storage_hold_reason(dict(status='pending_resume', storage_resume_requested_unix=99), 'first', record))
        self.assertIsNone(campaign.storage_hold_reason(dict(status='pending_resume', storage_resume_requested_unix=101), 'first', record))
        self.assertIsNone(campaign.storage_hold_reason(dict(status='pending'), 'other', record))

    def test_persistent_pressure_stops_owned_native_and_normal_exit_holds_restart(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); lane=root/'lane';lane.mkdir();out=lane/'result'
            runtime=root/'runtime';runtime.mkdir()
            solver=root/'solver'
            solver.write_text('''#!/usr/bin/python3
import json,signal,sys,time
from pathlib import Path
p=Path(sys.argv[1]);p.mkdir(exist_ok=True)
f=p/'starts';f.write_text(str(int(f.read_text())+1) if f.exists() else '1')
def stop(sig,frame):
 (p/'stopped').write_text('checkpoint preserved')
 raise SystemExit(0)
signal.signal(signal.SIGTERM,stop)
(p/'live').write_text('yes')
if (p/'allow_finish').exists():
 (p/'fullauto_status.json').write_text(json.dumps({'status':'converged'}))
else:
 while True:time.sleep(.01)
''');solver.chmod(0o755)
            manifest=lane/'manifest.json'
            manifest.write_text(json.dumps(dict(lane='test',gpus=['0'],environment={},binary=str(solver),
                binary_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(),poll_seconds=.05,
                jobs=[dict(id='one',output=str(out),arguments=[str(out)]),
                      dict(id='later',output=str(lane/'later'),arguments=[str(lane/'later')])])))
            wrapper=root/'run.py'
            wrapper.write_text('''import sys
from pathlib import Path
sys.path.insert(0,%r)
import run_fullauto_campaign as runner
import fullauto_storage as guard
runner.STORAGE_CHECK_INTERVAL_SECONDS=0
guard.BTRFS_PRESSURE_RECHECK_SECONDS=0
def health(lane,output=None):
 p=Path(output) if output else Path(lane)/'result'
 if (p/'live').exists() and not (p/'allow_finish').exists():
  return dict(ready=False,reason='persistent test pressure',btrfs=dict(ready=False,failure_kind='metadata_pressure',metadata_headroom_bytes=-458752))
 return dict(ready=True,reason='test ready')
guard.check_job_storage=health
try:runner.main()
except guard.StorageBlocked as e:print('HOLD storage:',e)
''' % str(SCRIPTS))
            environment=dict(os.environ,XDG_RUNTIME_DIR=str(runtime))
            def run():return subprocess.run([sys.executable,str(wrapper),str(manifest)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                            universal_newlines=True,env=environment,timeout=15)
            unrelated=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'],start_new_session=True)
            try:
                first=run();self.assertEqual(first.returncode,0,first.stderr)
                self.assertEqual((out/'starts').read_text(),'1');self.assertTrue((out/'stopped').exists())
                state=json.loads((lane/'campaign_status.json').read_text())
                self.assertEqual(state['jobs']['one']['status'],'storage_blocked')
                self.assertEqual(state['jobs']['one']['storage']['storage_confirmation']['attempts'],4)
                self.assertEqual(state['jobs']['one']['storage_stop']['action'],'terminated owned calculation group')
                self.assertNotIn('later',state['jobs']);self.assertIsNone(unrelated.poll())
                second=run();self.assertEqual(second.returncode,0,second.stderr)
                self.assertIn('HOLD one:',second.stdout);self.assertEqual((out/'starts').read_text(),'1')
                # A stale pending_resume plus an unrelated state-file update cannot clear the RAM hold.
                state['jobs']['one']['status']='pending_resume'
                state['unrelated_update']='new state-file mtime is not authorization'
                (lane/'campaign_status.json').write_text(json.dumps(state))
                held=run();self.assertEqual(held.returncode,0,held.stderr)
                self.assertIn('HOLD one:',held.stdout);self.assertEqual((out/'starts').read_text(),'1')
                # Explicit repair/unblock; preserve a previously completed other task.
                state['jobs']['one'].update(status='pending_resume',storage_resume_requested_unix=time.time()+1)
                state['jobs']['later']={'status':'converged','proof':'unchanged other task'}
                (lane/'campaign_status.json').write_text(json.dumps(state));(out/'allow_finish').write_text('yes')
                third=run();self.assertEqual(third.returncode,0,third.stderr)
                final=json.loads((lane/'campaign_status.json').read_text())
                self.assertEqual((out/'starts').read_text(),'2');self.assertEqual(final['jobs']['one']['status'],'converged')
                self.assertEqual(final['jobs']['later'],state['jobs']['later']);self.assertIsNone(unrelated.poll())
            finally:
                unrelated.terminate();unrelated.wait(timeout=5)


if __name__=='__main__':unittest.main()
