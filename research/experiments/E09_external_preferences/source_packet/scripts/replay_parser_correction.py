"""Re-evaluate frozen arrays under one text precision contract; no inference."""
from pathlib import Path
import subprocess,sys,json,time
R=Path(__file__).resolve().parents[1];attempts=[]
for name in ['evaluate_controlled.py','cross_assay.py','mapping_sensitivity.py',
             'evaluate_nogru_external.py','evaluate_hnr_component.py',
             'frequency_comparisons.py']:
    began=time.perf_counter()
    p=subprocess.run([sys.executable,'-X','utf8',str(R/'scripts'/name)],capture_output=True,text=True,encoding='utf-8')
    (R/'results'/('PARSER_REPLAY_'+name+'.log')).write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
    attempts.append(dict(script=name,exit_code=p.returncode,seconds=time.perf_counter()-began,no_model_inference=True))
    (R/'results'/'PARSER_REPLAY_EXECUTION.json').write_text(json.dumps(attempts,indent=2),encoding='utf-8')
    print(name,p.returncode,round(attempts[-1]['seconds'],2),flush=True)
    if p.returncode:raise SystemExit(p.returncode)
