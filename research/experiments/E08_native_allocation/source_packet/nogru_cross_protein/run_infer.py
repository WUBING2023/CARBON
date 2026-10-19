from pathlib import Path
import subprocess,shlex,json,time,hashlib,sys
ROOT=Path(__file__).resolve().parent
REMOTE='/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_cross_protein'
PY='/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3'
protocol=(ROOT/'PROTOCOL.md').read_bytes();code=(ROOT/'remote_infer.py').read_bytes()
import base64
payload=json.dumps({'source_b64':base64.b64encode(code).decode(),'protocol_b64':base64.b64encode(protocol).decode()}).encode()
boot='import sys,json,base64,pathlib,subprocess,os; r=pathlib.Path('+repr(REMOTE)+'); r.mkdir(parents=True,exist_ok=True); d=json.load(sys.stdin); (r/"PROTOCOL_TRAINING_FROZEN.md").write_bytes((r/"PROTOCOL.md").read_bytes()) if not (r/"PROTOCOL_TRAINING_FROZEN.md").exists() else None; (r/"remote_infer.py").write_bytes(base64.b64decode(d["source_b64"])); (r/"PROTOCOL.md").write_bytes(base64.b64decode(d["protocol_b64"])); os.environ["PYTHONDONTWRITEBYTECODE"]="1"; sys.exit(subprocess.call(['+repr(PY)+',str(r/"remote_infer.py")]))'
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join([PY,'-c',boot])]
t=time.perf_counter();p=subprocess.run(cmd,input=payload,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
(ROOT/'INFERENCE_LOG.jsonl').write_bytes(p.stdout)
(ROOT/'TRANSPORT_INFERENCE.json').write_text(json.dumps({'exit_code':p.returncode,'elapsed_s':time.perf_counter()-t,'code_sha256':hashlib.sha256(code).hexdigest(),'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'uploaded':'Only executable Python code and protocol; no local derived sequence inputs uploaded'},indent=2),encoding='utf-8')
print(p.stdout.decode('utf-8','replace'),flush=True)
sys.exit(p.returncode)
