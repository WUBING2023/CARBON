from pathlib import Path
import json,base64,subprocess,shlex
R=Path(__file__).resolve().parent;remote='/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr4';py='/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3';payload=json.dumps({'remote_nogru.py':base64.b64encode((R/'remote_nogru.py').read_bytes()).decode()}).encode()
boot="import sys,json,base64,pathlib,subprocess,os; r=pathlib.Path("+repr(remote)+"); d=json.load(sys.stdin); [(r/k).write_bytes(base64.b64decode(v)) for k,v in d.items()]; os.environ['PYTHONDONTWRITEBYTECODE']='1'; os.environ['OMP_NUM_THREADS']='4'; os.environ['OPENBLAS_NUM_THREADS']='4'; sys.exit(subprocess.call(["+repr(py)+",str(r/'remote_nogru.py')]))";cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join([py,'-c',boot])];p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=(R/'NOGRU_PROGRESS.txt').open('wb'));out,_=p.communicate(payload)
if p.returncode:raise RuntimeError((R/'NOGRU_PROGRESS.txt').read_text('utf8'))
d=json.loads(out)
for k,v in d['files'].items():(R/k).write_bytes(base64.b64decode(v))
print(json.dumps({k:v for k,v in d.items() if k!='files'}))
