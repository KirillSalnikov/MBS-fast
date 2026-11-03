"""A real subprocess verifies serial replay and terminal-job preservation."""
from pathlib import Path
import tempfile,json,hashlib,subprocess,unittest
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/run_fullauto_campaign.py'
class Resume(unittest.TestCase):
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
