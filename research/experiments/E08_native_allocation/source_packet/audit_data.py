from pathlib import Path
import subprocess,shlex,json,time
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import json,h5py,numpy as np
r=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');m=json.loads((r/'manifest.json').read_text())
out={'tasks':[],'previous_entries':[]}
for task,spec in m['tasks'].items():
 d=Path('/data/CARBON_PROJECT/数据/h5data')/task
 p=d/'_5ch'/(task+'_chr1_5ch.h5');raw=d/(task+'_chr1.h5')
 rec={'task':task,'counts':spec['counts'],'chr1_5ch_path':str(p),'chr1_raw_path':str(raw),'chr1_5ch_exists':p.exists(),'chr1_raw_exists':raw.exists()}
 if p.exists():
  with h5py.File(p,'r') as h:rec['h5_keys']={k:list(h[k].shape) for k in h.keys()};rec['x_finite']=bool(np.isfinite(h['X'][:256]).all())
  rec['bytes']=p.stat().st_size
 out['tasks'].append(rec)
for folder in ['functional_validation_20261002','expanded_functional_validation_20261003','independent_context_validation_20261003','interpretability_strengthening_20261003','functional_candidates_20260930','可解释性分析','比较实验']:
 d=Path('/data/CARBON_PROJECT')/folder
 out['previous_entries'].append({'directory':str(d),'files':[str(x.relative_to(d)) for x in d.iterdir()] if d.exists() else []})
print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True)
if p.returncode: raise RuntimeError(p.stderr.decode('utf-8','replace'))
(ROOT/'DATA_AUDIT.json').write_bytes(p.stdout);print(p.stdout.decode('utf-8'));print('Elapsed seconds',time.perf_counter()-t)
