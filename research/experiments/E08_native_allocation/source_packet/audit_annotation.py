from pathlib import Path
import subprocess,shlex,json,time
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import json,os
roots=['/data/CARBON_PROJECT/数据','/data/CARBON_PROJECT/数据处理','/data/CARBON_PROJECT/independent_context_validation_20261003','/data/CARBON_PROJECT/interpretability_strengthening_20261003','/data/CARBON_PROJECT/expanded_functional_validation_20261003']
out={}
for root in roots:
 files=[]
 for base,dirs,names in os.walk(root):
  depth=len(Path(base).relative_to(root).parts)
  if depth>=4:dirs.clear()
  for name in names:
   if root.endswith(('数据','数据处理')):
    keep=name.endswith(('.gtf','.gtf.gz','.gff','.gff3','.bed','.fa.fai','.bigWig','.bw')) and ('gencode' in name.lower() or 'phylo' in name.lower() or name.endswith('.gtf'))
   else:keep=name.endswith(('.py','.json','.csv','.tsv')) and any(k in name.lower() for k in ('input','manifest','cohort','case','coord','meta','run','predict'))
   if keep:
    p=Path(base)/name;files.append({'path':str(p),'bytes':p.stat().st_size})
 out[root]=files
print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
p=subprocess.run(cmd,capture_output=True,check=True);(ROOT/'ANNOTATION_AUDIT.json').write_bytes(p.stdout);print(p.stdout.decode('utf-8'))
