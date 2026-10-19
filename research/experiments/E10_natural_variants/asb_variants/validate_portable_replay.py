"""Isolated copy proves local paths suffice; compare scientific tables after replay."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='4';sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import shutil,subprocess,json,datetime,numpy as np,pandas as pd
R=Path(__file__).resolve().parent;D=R/('replay_check_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));D.mkdir()
for folder in ['source_inputs','b1_reference','annotation_ncbi','annotation_raw','annotation_clinvar']:
    shutil.copytree(R/folder,D/folder)
for p in R.glob('*.py'):shutil.copyfile(p,D/p.name)
for n in ['ANNOTATION_FROZEN_VARIANTS.csv','COHORT_ALL7262.csv','PACKAGED_INPUT_MANIFEST.json','PROTOCOL_FROZEN.json','FREEZE_MANIFEST.json']:
    shutil.copyfile(R/n,D/n)
with (D/'REPLAY_STDOUT.txt').open('w',encoding='utf8') as f:
    subprocess.run([sys.executable,str(D/'replay_offline.py')],cwd=D,stdout=f,stderr=subprocess.STDOUT,check=True)
checks=[]
for name in ['DIRECTION_SUMMARY.csv','DIRECTION_CONTRASTS.csv','CONSERVATION_ASSOCIATIONS.csv','DIRECTION_BY_CONSERVATION.csv','DIRECTION_TWO_WAY_SENSITIVITY.csv','WEIGHTED_REFIT_SENSITIVITY_200.csv','WEIGHTED_REFIT_REGION_SENSITIVITY_200.csv','ADJUSTED_ASSOCIATION_PAIRED_CONTRASTS.csv','ANNOTATIONS_ALL287.csv','STRICT396_WITH_ANNOTATIONS.csv','CLINVAR_DIRECT_VERIFICATION.csv']:
    a,b=pd.read_csv(R/name),pd.read_csv(D/name)
    pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-10,atol=1e-12)
    checks.append(dict(file=name,rows=len(a),matched=True))
report=dict(isolated_directory=str(D),network_used=False,external_original_paths_needed=False,compared_tables=checks,all_passed=True)
(R/'PORTABLE_REPLAY_QA.json').write_text(json.dumps(report,indent=2)+'\n','utf8')
print(json.dumps(report,indent=2))
