from pathlib import Path
import subprocess,json,hashlib,datetime,time
R=Path(__file__).resolve().parent
payload=(R/'MAPPING_PAYLOAD.json').read_text('utf8');remote=(R/'remote_map_fr2.py').read_text('utf8')
program='import sys,io,json\n_payload='+repr(payload)+'\nsys.stdin=io.StringIO(_payload)\nexec('+repr(remote)+')\n'
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15','PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 /usr/bin/python3 -']
start=time.perf_counter();r=subprocess.run(cmd,input=program.encode('utf8'),capture_output=True,timeout=300)
(R/'MAPPING_REMOTE_LOG.txt').write_text(r.stderr.decode('utf8',errors='replace'),encoding='utf8')
if r.returncode:raise RuntimeError(r.stderr.decode('utf8',errors='replace'))
data=json.loads(r.stdout);data['local_provenance']={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CPU_wall_seconds':time.perf_counter()-start,'payload_sha256':hashlib.sha256(payload.encode()).hexdigest(),'script_sha256':hashlib.sha256(remote.encode()).hexdigest(),'model_queries':0,'readonly_public_reference_mapping':True}
(R/'REFERENCE_MAPPING.json').write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps({'source_rows':len(data['records']),'unique_map':sum(r['unique_natural_human_RNA_mapping'] for r in data['records']),'complete_Phi':sum(r['unique_natural_human_RNA_mapping'] and r['maps'][0]['phi_missing']==0 for r in data['records']),'CPU_wall_seconds':data['local_provenance']['CPU_wall_seconds']}))
