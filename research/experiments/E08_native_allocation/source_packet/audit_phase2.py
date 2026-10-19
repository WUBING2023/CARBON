from pathlib import Path
import subprocess,shlex,json,time
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import json,h5py
r=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');m=json.loads((r/'manifest.json').read_text())
out={'tasks':[],'source_files':{},'no_gru_source_files':{}}
for task in m['tasks']:
 d=Path('/data/CARBON_PROJECT/数据/h5data')/task
 for ch in ['X','Y','M']:
  p=d/'_5ch'/(task+'_chr'+ch+'_5ch.h5');raw=d/(task+'_chr'+ch+'.h5')
  rec={'task':task,'chromosome':ch,'path':str(p),'raw_path':str(raw),'exists':p.exists(),'raw_exists':raw.exists()}
  if p.exists():
   with h5py.File(p,'r') as h:rec['shapes']={k:list(h[k].shape) for k in h.keys()}
  out['tasks'].append(rec)
for name in ['matched_train.py','matched_common.py','source/MARCP/models/marcp_model.py']:
 out['source_files'][name]=(r/name).read_text()
d=Path('/data/CARBON_PROJECT/消融实验/全量消融运行/workspace/MARCP')
for p in d.rglob('*.py'):
 if 'no_gru' in p.name.lower() or (p.name=='marcp_model.py'):
  out['no_gru_source_files'][str(p)]=p.read_text()
print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True)
if p.returncode:raise RuntimeError(p.stderr.decode('utf-8','replace'))
(ROOT/'PHASE2_AUDIT.json').write_bytes(p.stdout)
j=json.loads(p.stdout);print(json.dumps({'tasks':j['tasks'],'no_gru_source_files':list(j['no_gru_source_files'])},ensure_ascii=False,indent=2));print('seconds',time.perf_counter()-t)
