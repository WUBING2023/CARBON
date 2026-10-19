from pathlib import Path
import json,base64,subprocess,shlex,time
R=Path(__file__).resolve().parent;remote='/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr4';py='/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3'
files={'WINDOWS.csv':R/'WINDOWS.csv','remote_fr4.py':R/'remote_fr4.py','remote_extract.py':R/'remote_extract.py','PROTOCOL_FR4.md':R.parent/'PROTOCOL_FR4.md','PROTOCOL_FR4_HOLDOUT_ADDENDUM.md':R.parent/'PROTOCOL_FR4_HOLDOUT_ADDENDUM.md'}
payload=json.dumps({k:base64.b64encode(p.read_bytes()).decode() for k,p in files.items()}).encode()
boot="import sys,json,base64,pathlib,subprocess,os; r=pathlib.Path("+repr(remote)+"); r.mkdir(parents=True,exist_ok=True); d=json.load(sys.stdin); [(r/k).write_bytes(base64.b64decode(v)) for k,v in d.items()]; os.environ['PYTHONDONTWRITEBYTECODE']='1'; os.environ['OMP_NUM_THREADS']='4'; os.environ['OPENBLAS_NUM_THREADS']='4'; sys.exit(subprocess.call(["+repr(py)+",str(r/'remote_fr4.py')]))"
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join([py,'-c',boot])]
t=time.perf_counter();p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=(R/'REMOTE_PROGRESS.jsonl').open('wb'));out,_=p.communicate(payload)
if p.returncode:raise RuntimeError((R/'REMOTE_PROGRESS.jsonl').read_text('utf8'))
d=json.loads(out)
for k,v in d['files'].items():(R/k).write_bytes(base64.b64decode(v))
print(json.dumps({'local_wall_seconds':time.perf_counter()-t,**{k:d['summary'][k] for k in ['wall_seconds','queries','native_valid_windows','native_excluded_events']}}))
