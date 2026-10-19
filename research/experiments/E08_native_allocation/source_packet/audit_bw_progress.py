from pathlib import Path
import subprocess,shlex,json
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import json,os
O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_matched')
out={'training_manifest':json.loads((O/'TRAINING_MANIFEST.json').read_text()) if (O/'TRAINING_MANIFEST.json').exists() else None,'progress':[],'data_directories':[],'bw_paths':[]}
for p in (O/'outputs').glob('*/*'):
 rec={'directory':str(p)}
 for f in ['progress.json','process.json','complete.json']:
  if (p/f).exists():rec[f]=json.loads((p/f).read_text())
 if (p/'console.log').exists():rec['console_tail']=(p/'console.log').read_text().splitlines()[-2:]
 out['progress'].append(rec)
d=Path('/data/CARBON_PROJECT/数据')
out['data_directories']=[str(p) for p in d.iterdir() if p.is_dir()]
for start in [d,Path('/data/CARBON_PROJECT/数据处理')]:
 for root,dirs,files in os.walk(start):
  if len(Path(root).relative_to(start).parts)>=4:dirs[:]=[]
  for f in files:
   if f.lower().endswith(('.bw','.bigwig')):out['bw_paths'].append(str(Path(root)/f))
print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
p=subprocess.run(cmd,capture_output=True)
if p.returncode:raise RuntimeError(p.stderr.decode())
(ROOT/'BW_TRAINING_AUDIT.json').write_bytes(p.stdout);j=json.loads(p.stdout);print(json.dumps({'progress':j['progress'],'data_directories':j['data_directories'],'bw_count':len(j['bw_paths']),'bw_example':j['bw_paths'][:10]},ensure_ascii=False,indent=2))
