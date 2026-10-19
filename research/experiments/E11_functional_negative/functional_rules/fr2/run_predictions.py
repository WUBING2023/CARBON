from pathlib import Path
import subprocess,json,base64,shlex,datetime,time
R=Path(__file__).resolve().parent;remote='/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr2';py='/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3'
files={'INPUTS.npz':R/'INPUTS.npz','PROTOCOL_FR2.md':R.parent/'PROTOCOL_FR2.md','remote_predict_fr2.py':R/'remote_predict_fr2.py'}
payload=json.dumps({k:base64.b64encode(p.read_bytes()).decode() for k,p in files.items()}).encode()
boot="import sys,json,base64,pathlib,subprocess,os; r=pathlib.Path("+repr(remote)+"); r.mkdir(parents=True,exist_ok=True); d=json.load(sys.stdin); [(r/k).write_bytes(base64.b64decode(v)) for k,v in d.items()]; os.environ['PYTHONDONTWRITEBYTECODE']='1'; os.environ['OMP_NUM_THREADS']='4'; os.environ['OPENBLAS_NUM_THREADS']='4'; sys.exit(subprocess.call(["+repr(py)+",str(r/'remote_predict_fr2.py')]))"
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join([py,'-c',boot])]
start=time.perf_counter();p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=(R/'REMOTE_PROGRESS.jsonl').open('wb'));out,_=p.communicate(payload)
if p.returncode:raise RuntimeError((R/'REMOTE_PROGRESS.jsonl').read_text('utf8'))
d=json.loads(out);data=base64.b64decode(d['arrays_base64']);(R/'PREDICTIONS.npz').write_bytes(data);(R/'RESULT.json').write_text(json.dumps(d['result'],indent=2),encoding='utf8')
print(json.dumps({'exit_code':p.returncode,'local_wall_seconds':time.perf_counter()-start,**{k:d['result'][k] for k in ['model_queries','n_rows','wall_seconds','predictions_sha256']}}))
