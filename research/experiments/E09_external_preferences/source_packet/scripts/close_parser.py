"""Snapshot and apply a uniform literal numeric text parsing contract."""
from pathlib import Path
import shutil, hashlib, json
R=Path(__file__).resolve().parents[1]
old=R/'revisions'/'parser_default_20261004'
old.mkdir(parents=True,exist_ok=True)
names=['evaluate_rules.py','cross_assay.py','mapping_sensitivity.py','evaluate_wide.py',
       'evaluate_nogru_external.py','evaluate_hnr_component.py','hnrnpk_stability.py',
       'frequency_comparisons.py','plot_external.py','package_packet.py']
records=[]
for folder,selected in [('scripts',[R/'scripts'/n for n in names]),
                         ('results',[p for p in (R/'results').glob('*') if p.is_file() and p.suffix in ['.csv','.json']]),
                         ('root',[R/'MANIFEST.json',R/'QA.json'])]:
    d=old/folder;d.mkdir(exist_ok=True)
    for p in selected:
        dest=d/p.name
        if not dest.exists():shutil.copy2(p,dest)
        records.append(dict(path=f'{folder}/{p.name}',sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
(old/'MANIFEST.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
for name in names:
    p=R/'scripts'/name;s=p.read_text(encoding='utf-8')
    if 'pd.read_csv(' in s:
        s=s.replace('pd.read_csv(', 'read_csv_literal(')
        s='from literal_tables import read_csv_literal\n'+s
        p.write_text(s,encoding='utf-8')
print(json.dumps(dict(snapshot_files=len(records),uniform_literal_readers=names),indent=2))
