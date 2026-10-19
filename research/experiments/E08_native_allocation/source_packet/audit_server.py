from pathlib import Path
import subprocess,shlex,json,time
ROOT=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import os,json,subprocess,platform
r=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
m=json.loads((r/'manifest.json').read_text())
out={'hostname':platform.node(),'cpus':os.cpu_count(),'loadavg':os.getloadavg(),
 'gpu':subprocess.check_output(['nvidia-smi','--query-gpu=index,name,memory.used,memory.total,utilization.gpu','--format=csv,noheader'],text=True),
 'gpu_jobs':subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_gpu_memory','--format=csv,noheader'],text=True),
 'manifest':m,'matched_files':[str(x.relative_to(r)) for x in r.iterdir()],
 'rbp_entries':[{'name':x.name,'directory':x.is_dir()} for x in Path('/data/CARBON_PROJECT').iterdir()],
 'checkpoint_files':[{'path':str(x),'bytes':x.stat().st_size} for ext in ['*.pt','*.pth'] for x in r.rglob(ext)],
 'py_files':[str(x) for x in r.rglob('*.py')],
 'memory':Path('/proc/meminfo').read_text().splitlines()[:5]}
print(json.dumps(out,ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True)
if p.returncode: raise RuntimeError(p.stderr.decode('utf-8','replace'))
(ROOT/'SERVER_AUDIT.json').write_bytes(p.stdout)
o=json.loads(p.stdout)
print(json.dumps({'elapsed_s':time.perf_counter()-t,'cpus':o['cpus'],'loadavg':o['loadavg'],'gpu':o['gpu'],'gpu_jobs':o['gpu_jobs'],'tasks':list(o['manifest']['tasks']),'checkpoints':len(o['checkpoint_files']),'matched_files':o['matched_files'],'py_files':o['py_files'],'rbp_dirs':[x['name'] for x in o['rbp_entries'] if x['directory']]},ensure_ascii=False,indent=2))
