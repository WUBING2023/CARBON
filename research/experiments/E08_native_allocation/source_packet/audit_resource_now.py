from pathlib import Path
import subprocess,shlex,json
ROOT=Path(__file__).resolve().parent
CODE=r'''
import subprocess,json,time
q=lambda args:subprocess.run(args,capture_output=True,text=True).stdout
gpu=q(['nvidia-smi','--query-gpu=index,uuid,name,memory.used,memory.total,utilization.gpu','--format=csv,noheader'])
jobs=q(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,process_name,used_gpu_memory','--format=csv,noheader'])
pids=[v.split(',')[1].strip() for v in jobs.splitlines() if ',' in v]
print(json.dumps({'unix':time.time(),'gpu':gpu,'jobs':jobs,'processes':{p:q(['ps','-p',p,'-o','user=,pid=,etime=,args=']) for p in pids},'load':q(['cat','/proc/loadavg'])},ensure_ascii=False,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
p=subprocess.run(cmd,capture_output=True);print(p.stdout.decode());(ROOT/'RESOURCE_1143.json').write_bytes(p.stdout)
if p.returncode:raise RuntimeError(p.stderr.decode())
