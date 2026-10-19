from pathlib import Path
import subprocess,json,shlex
R=Path(__file__).resolve().parent
code="""import os,json,subprocess,datetime
root='/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/'
own=[]
for name in os.listdir('/proc'):
 if not name.isdigit() or int(name)==os.getpid():continue
 try:
  args=open('/proc/'+name+'/cmdline','rb').read().decode(errors='replace').split('\\0')
  if len(args)>1 and 'python' in os.path.basename(args[0]) and args[1].startswith(root):own.append({'pid':int(name),'script':args[1]})
 except (OSError,UnicodeError):pass
gpu=subprocess.run(['nvidia-smi','--query-gpu=index,uuid,memory.used,utilization.gpu','--format=csv,noheader'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
g2=[v.strip() for v in gpu.stdout.splitlines() if v.split(',')[0].strip()=='2']
print(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'own_active_python_processes':own,'gpu2_snapshot':g2,'gpu_query_exit':gpu.returncode,'own_batch_finished_utc':'2026-10-04T13:34:10+00:00','no_resource_reservation_remaining':not own}))
"""
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['python3','-c',code])]
p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=45)
if p.returncode:raise RuntimeError('Read-only SSH resource check failed: '+p.stderr.decode(errors='replace'))
d=json.loads(p.stdout);(R/'RESOURCE_RELEASE.json').write_text(json.dumps(d,indent=2),encoding='utf-8');print(json.dumps(d))
