from pathlib import Path
import subprocess,shlex,json
R=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import os,json,time
out={'bw':[],'task_directories':[],'start_unix':time.time()}
for start in [Path('/data/CARBON_PROJECT/数据'),Path('/data/CARBON_PROJECT/数据处理')]:
 for root,dirs,files in os.walk(start):
  if len(Path(root).relative_to(start).parts)>=8:dirs[:]=[]
  if 'PTBP1' in root and len(out['task_directories'])<100:out['task_directories'].append(root)
  for f in files:
   if f.lower().endswith(('.bw','.bigwig')):out['bw'].append(str(Path(root)/f))
out['seconds']=time.time()-out['start_unix'];print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
p=subprocess.run(cmd,capture_output=True)
if p.returncode:raise RuntimeError(p.stderr.decode())
(R/'PTBP_BW_SEARCH.json').write_bytes(p.stdout);j=json.loads(p.stdout);print(json.dumps({'ptbp_bw':[v for v in j['bw'] if 'PTBP' in v],'bw_files':len(j['bw']),'seconds':j['seconds'],'directories':j['task_directories'][:25]},ensure_ascii=False,indent=2))
