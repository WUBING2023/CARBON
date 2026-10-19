#!/usr/bin/env python3
"""Evaluate the included AARS checkpoint on CPU and compare six saved demo metrics.
This checks the included demonstration, not historical MARCP equivalence or training.
"""
import argparse,csv,json,subprocess,sys,hashlib,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for stream in (sys.stdout,sys.stderr):
 if hasattr(stream,'reconfigure'):stream.reconfigure(encoding='utf-8',errors='replace')
KEYS=['Mean_PCC','Spearman','Signal_AUC','Signal_AP','Peak_AUC','Peak_AP']
def load(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return {r['Metric']:float(r['Value']) for r in csv.DictReader(f,delimiter='\t')}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'outputs/demo_cpu');p.add_argument('--tolerance',type=float,default=1e-5);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
 checkpoint=ROOT/'carbon_trained_model/AARS_K562_ENCSR825SVO_demo/model_params.pth';data=ROOT/'data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr3_5ch.h5'
 args=['carbon_train.py','--eval_only','--load_model',str(checkpoint),'--CUDA','-1','--signalMax','0','--seqLength','101','--bsz','64','--encode','OneHot','--modelFramework','Treat','--input_channels','5','-test',str(data),'--task','AARS_K562_ENCSR825SVO','--lossM','CombinedProfileLoss','--alpha','1.42','--beta','0.54','--gamma','1.80','--output',str(out)]
 child_env=dict(os.environ,PYTHONUTF8="1",PYTHONIOENCODING="utf-8")
 result=subprocess.run([sys.executable,"-X","utf8",*args],cwd=ROOT,env=child_env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
 (out/'evaluation.log').write_text(result.stdout,encoding='utf-8');print(result.stdout[-3000:])
 if result.returncode:raise SystemExit(result.returncode)
 ref=load(ROOT/'carbon_trained_model/AARS_K562_ENCSR825SVO_demo/test_metrics_summary.txt');actual=load(out/'test_metrics_summary.txt');diff={k:abs(ref[k]-actual[k]) for k in KEYS};maximum=max(diff.values())
 report=dict(scope='Included AARS checkpoint CPU evaluation; no training or MARCP equivalence proof.',test_windows=1283,tolerance=a.tolerance,maximum_absolute_error=maximum,passed=maximum<=a.tolerance,reference={k:ref[k] for k in KEYS},actual={k:actual[k] for k in KEYS},absolute_errors=diff,checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),test_h5_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),python=sys.version)
 (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)
if __name__=='__main__':main()
