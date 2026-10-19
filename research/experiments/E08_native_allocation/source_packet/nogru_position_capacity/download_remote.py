from pathlib import Path
import subprocess,shlex,json,hashlib,time,io,zipfile
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import io,zipfile,sys
r=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_position_capacity')
buf=io.BytesIO()
with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as z:
 for p in r.rglob('*'):
  if p.is_file() and (p.suffix in ('.json','.jsonl','.npz','.py','.md') or p.name=='best.pth'):z.write(p,p.relative_to(r))
sys.stdout.buffer.write(buf.getvalue())
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True,check=True)
archive=ROOT/'remote_results.zip';archive.write_bytes(p.stdout)
with zipfile.ZipFile(io.BytesIO(p.stdout)) as z:
 assert all(not n.startswith('/') and '..' not in Path(n).parts for n in z.namelist())
 z.extractall(ROOT/'remote')
(ROOT/'TRANSPORT_DOWNLOAD.json').write_text(json.dumps({'bytes':len(p.stdout),'sha256':hashlib.sha256(p.stdout).hexdigest(),'elapsed_s':time.perf_counter()-t,'source':'Task-owned new server directory only','stderr':p.stderr.decode('utf-8','replace')},indent=2),encoding='utf-8')
print(json.dumps({'bytes':len(p.stdout),'seconds':time.perf_counter()-t,'archive':str(archive)}))
