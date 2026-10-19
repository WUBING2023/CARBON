#!/usr/bin/env python3
"""Verify packaged file identities, alias targets and current manuscript bindings."""
from pathlib import Path
import argparse,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'outputs/package_verification.json');a=p.parse_args();failures=[];manifest=rows(ROOT/'research/catalog/SOURCE_MANIFEST.csv')
 for r in manifest:
  q=ROOT/r['path']
  if not q.is_file() or sha(q)!=r['sha256']:failures.append(r['path'])
  approved=(r['path']=='paper/source_data/S23_external_order_preferences/original/sources/Zscores_420_origrncmpt.tar.gz' and int(r['bytes'])==51769595 and r['sha256']=='bf65e4bdd12fb1f8fd5db585e5165db24cc1dba6a1b2708f2cc0e12255a676bd')
  if int(r['bytes'])>50000000 and not approved:failures.append('Oversized candidate: '+r['path'])
 aliases=rows(ROOT/'research/catalog/ASSET_ALIASES.csv')
 for r in aliases:
  q=ROOT/r['canonical']
  if not q.is_file() or sha(q)!=r['sha256']:failures.append('Alias target: '+r['path'])
 binding=json.loads((ROOT/'paper/MANUSCRIPT_BINDING.json').read_text())
 for r in binding['files']:
  if sha(ROOT/next((x['canonical'] for x in aliases if x['path']==r['path']),r['path']))!=r['sha256']:failures.append('Manuscript binding: '+r['path'])
 for r in rows(ROOT/'paper/DATASET_DECLARATIONS.csv'):
  if r['exists']!='True':failures.append('Declared data directory: '+r['declared_path'])
 datasets=rows(ROOT/'paper/SUPPLEMENT_DATA_BINDINGS.csv')
 if sorted(r['dataset'] for r in datasets)!=sorted('S'+str(n) for n in range(1,28)):failures.append('Expected each dataset S1-S27 exactly once')
 amap={r['path']:r['canonical'] for r in aliases}
 choices=json.loads((ROOT/'paper/PRIMARY_DATA_SELECTIONS.json').read_text(encoding='utf-8'))['selections'];policy={r['dataset']:r for r in choices}
 if len(choices)!=27 or set(policy)!={'S'+str(n) for n in range(1,28)}:failures.append('Expected explicit scientific selections S1-S27')
 for r in datasets:
  if r['status']!='bound':failures.append('Supplementary dataset: '+r['dataset'])
  n=int(r['dataset'][1:]);prefix='paper/source_data/S1_S14_original/' if n<=14 else 'paper/source_data/S'+str(n)+'_'
  declared=r.get('declared_path','').rstrip('/');selected=r.get('declared_primary_path','')
  if not selected.startswith(prefix) or not (selected==declared or selected.startswith(declared+'/')):failures.append('Dataset declared boundary mismatch: '+r['dataset'])
  choice=policy.get(r['dataset'],{});actual=ROOT/r['primary_path'];source=ROOT/amap.get(selected,selected)
  if selected!=choice.get('declared_primary_path') or r['primary_path']!=choice.get('primary_path') or r.get('primary_role')!=choice.get('primary_role'):failures.append('Scientific primary selection mismatch: '+r['dataset'])
  if not source.is_file():failures.append('Declared primary unavailable: '+r['dataset'])
  elif actual.is_file() and sha(source)!=sha(actual):failures.append('Declared/canonical content mismatch: '+r['dataset'])
  if not actual.is_file():failures.append('Dataset primary file absent: '+r['dataset']+' '+r['primary_path'])
  elif sha(actual)!=r.get('primary_sha256'):failures.append('Dataset primary SHA mismatch: '+r['dataset'])
  if n==27 and not r['primary_path'].endswith('/ALL66_REPORTER_BIOLOGICAL_MEASUREMENTS.csv'):failures.append('Reporter primary mismatch')
 for r in rows(ROOT/'paper/FIGURE_BINDINGS.csv'):
  if not r['path'] or not (ROOT/r['path']).is_file() or sha(ROOT/r['path'])!=r['sha256']:failures.append('Figure output: '+r['reference'])
 result=dict(passed=not failures,checked_files=len(manifest),alias_targets=len(aliases),primary_data_files=len(datasets),bound_paper_files=len(binding['files']),failures=failures,scope='File identity, candidate size and binding integrity; no upstream reproducibility claim.')
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
if __name__=='__main__':main()
