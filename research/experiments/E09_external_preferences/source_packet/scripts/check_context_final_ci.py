"""Independent final contrast-name and multiplicity interval readback."""
from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];A=R.parent/'01_context'/'nogru_position_capacity'
z=np.load(A/'bootstrap'/'all_core_adjusted_PTBP1.npz');names=list(z['methods']);m=z['means']
saved=pd.read_csv(A/'CONTRASTS.csv',float_precision='round_trip')
saved=saved.query("scene=='all_core' and cohort=='adjusted' and group=='PTBP1' and seed=='ensemble'")
rows=[]
for row in saved.itertuples():
    if row.comparison=='full_minus_control':left='aligned_ensemble';right=row.control+'_ensemble'
    elif row.comparison=='position_capacity_minus_nogru':left='nogru_position_capacity_ensemble';right='nogru_ensemble'
    else:raise AssertionError(row.comparison)
    d=m[:,names.index(left)]-m[:,names.index(right)]
    for family in [1,2,3,15,45,90]:
        ci=np.quantile(d,[.025/family,1-.025/family]);expected=[getattr(row,f'bonf{family}_low'),getattr(row,f'bonf{family}_high')]
        assert np.max(np.abs(ci-expected))<1e-11,(row.comparison,row.control,family)
        rows.append(dict(comparison=row.comparison,control=row.control,interval_family=family,
                         independent_low=float(ci[0]),independent_high=float(ci[1]),
                         max_error=float(np.max(np.abs(ci-expected)))))
pd.DataFrame(rows).to_csv(R/'results'/'CROSS_AUDIT_01_PC_FINAL_CI.csv',index=False)
(R/'results'/'CROSS_AUDIT_01_PC_FINAL_CI_QA.json').write_text(json.dumps(dict(n_checks=len(rows),
    all_passed=True,explicit_comparison_field_used=True,bootstrap_source='actual paired gene saved arm arrays',
    original_10000_draws_already_independently_reproduced=True,max_abs_error=max(v['max_error'] for v in rows)),indent=2))
print(pd.DataFrame(rows).query("control=='nogru_position_capacity' and interval_family==2").to_string(index=False))
print(len(rows),'final contrast/interval checks passed')
