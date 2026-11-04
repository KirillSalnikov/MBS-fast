"""A real subprocess verifies serial replay and terminal-job preservation."""
from pathlib import Path
import tempfile,json,hashlib,subprocess,unittest,importlib.util,sys
from unittest.mock import patch
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/run_fullauto_campaign.py'
sys.path.insert(0,str(SCRIPT.parent))
spec=importlib.util.spec_from_file_location('campaign',SCRIPT)
campaign=importlib.util.module_from_spec(spec);spec.loader.exec_module(campaign)
class Resume(unittest.TestCase):
 def test_storage_pause_is_visible_without_discarding_resume_metadata(self):
  statuses={'job':{'status':'running','pid':123,'command':['solver'],'elapsed_seconds':42}}
  health={'ready':False,'reason':'inode reserve'}
  calls=[]
  with patch.object(campaign,'record_runtime') as runtime:
   campaign.publish_storage_block('/lane',statuses,'job',health,lambda *args:calls.append(args))
  self.assertEqual(calls,[('storage_blocked','job')])
  self.assertEqual(statuses['job']['status'],'storage_blocked')
  self.assertEqual(statuses['job']['elapsed_seconds'],42)
  self.assertEqual(statuses['job']['command'],['solver'])
  runtime.assert_called_once()
 def test_runtime_diagnostic_survives_queue_write_failure(self):
  def fail(*args):raise OSError(28,'full')
  with patch.object(campaign,'record_runtime') as runtime:
   campaign.publish_storage_block('/lane',{},'job',{'ready':False},fail)
  self.assertEqual(runtime.call_args.args[1]['status'],'storage_blocked')
 def test_terminal_and_interrupted(self):
  with tempfile.TemporaryDirectory()as d:
   p=Path(d);exe=p/'solver';out=p/'result'
   exe.write_text('#!/usr/bin/python3\nfrom pathlib import Path\nimport sys,json\np=Path(sys.argv[1]);p.mkdir(exist_ok=True)\nf=p/"calls";f.write_text(str(int(f.read_text())+1) if f.exists() else "1")\n(p/"fullauto_status.json").write_text(json.dumps({"status":"converged"}))\n');exe.chmod(0o755)
   m={'lane':'test','gpus':['0'],'environment':{},'binary':str(exe),'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'poll_seconds':.05,'jobs':[{'id':'job','output':str(out),'arguments':[str(out)]}]};manifest=p/'manifest.json';manifest.write_text(json.dumps(m))
   def run():subprocess.run(['/usr/bin/python3',str(SCRIPT),str(manifest)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10)
   run();s=json.loads((p/'campaign_status.json').read_text());self.assertEqual(s['jobs']['job']['status'],'converged');run();self.assertEqual((out/'calls').read_text(),'1')
   j=s['jobs']['job'];j.update(status='running',pid=99999999);j.pop('finished_unix',None);j.pop('returncode',None);(p/'campaign_status.json').write_text(json.dumps(s))
   run();s=json.loads((p/'campaign_status.json').read_text());self.assertEqual((out/'calls').read_text(),'2');self.assertEqual(s['jobs']['job']['status'],'converged');self.assertEqual(s['jobs']['job']['attempt'],2)
if __name__=='__main__':unittest.main()
