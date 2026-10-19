from pathlib import Path
import subprocess,shlex,json
O=Path(__file__).resolve().parent
CODE=r'''
from pathlib import Path
import json,h5py,numpy as np,pandas as pd,ast,time
P=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context');rec=[];start=time.time()
def signal_parse(value):
 if not isinstance(value,str):return value
 class NumericWrapper(ast.NodeTransformer):
  def visit_Call(self,node):
   ok=(isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id in ['np','numpy'] and node.func.attr in ['float32','float64','int32','int64'])
   assert ok and len(node.args)==1 and not node.keywords,'Unexpected constructor in numeric profile literal'
   return self.visit(node.args[0])
 tree=NumericWrapper().visit(ast.parse(value,mode='eval'));return ast.literal_eval(tree)
for root in [P,P/'phase2_ambiguity']:
 m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
 for s in m['tasks']:
  with h5py.File(s['input_path'],'r') as h:y=np.asarray(h['Y']);b=np.asarray(h['Y_ctl']);x=np.asarray(h['X'])
  assert np.isfinite(y).all() and np.isfinite(b).all() and np.isfinite(x).all() and (y>=0).all() and (b>=0).all()
  df=pd.read_hdf(s['raw_path'],key='df');df.columns=['chr','start','end','peak','score','strand','signal','input','seq'];df=df[~df['seq'].str.contains('N',case=False)].dropna();assert len(df)==len(y)
  ix=np.unique(np.linspace(0,max(0,len(df)-1),min(64,len(df)),dtype=int));diff=[];reverse=[]
  for i in ix:
   v=df.iloc[i];a=np.asarray(signal_parse(v['signal'])).ravel();c=np.asarray(signal_parse(v['input'])).ravel();assert len(a)==len(c)==101
   diff.append(max(float(np.max(abs(y[i]-a))),float(np.max(abs(b[i]-c)))));reverse.append(max(float(np.max(abs(y[i]-a[::-1]))),float(np.max(abs(b[i]-c[::-1])))))
  rec.append({'task':s['task'],'chromosome':s.get('chromosome','1'),'windows':len(y),'finite_nonnegative':True,'raw_profile_QA_samples':len(ix),'raw_direct_max_error':max(diff) if diff else None,'raw_reverse_max_error':max(reverse) if reverse else None,'raw_direct_all_exact':all(v==0 for v in diff),'native_unit':'Stored native profile units; not unverified rawread counts'})
out={'seconds':time.time()-start,'source_checks':rec};(P/'NATIVE_SOURCE_QA.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',CODE])]
p=subprocess.run(cmd,capture_output=True)
if p.returncode:raise RuntimeError(p.stderr.decode())
(O/'NATIVE_SOURCE_QA.json').write_bytes(p.stdout);print(p.stdout.decode())
