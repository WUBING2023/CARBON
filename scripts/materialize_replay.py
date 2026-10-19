#!/usr/bin/env python3
"""Copy this repository to new scratch and restore SHA-checked duplicate asset paths."""
from pathlib import Path
import argparse,csv,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--destination',required=True,type=Path);a=p.parse_args();dest=a.destination.resolve()
 if dest==ROOT.resolve() or ROOT.resolve() in dest.parents:raise SystemExit('Choose a new destination outside this repository.')
 if dest.exists():raise SystemExit('Destination already exists; choose a new scratch directory.')
 def ignore(src,names):
  skipped={n for n in names if n in {'.git','__pycache__'} or n.endswith('.pyc')}
  if Path(src).resolve()==ROOT.resolve():skipped.add('outputs')
  return skipped
 shutil.copytree(ROOT,dest,ignore=ignore)
 with (ROOT/'research/catalog/ASSET_ALIASES.csv').open(encoding='utf-8',newline='') as f:aliases=list(csv.DictReader(f))
 for r in aliases:
  src=dest/r['canonical'];target=dest/r['path'];target.resolve().relative_to(dest)
  if hashlib.sha256(src.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('Canonical asset checksum mismatch: '+r['canonical'])
  target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,target)
 print('Materialized',len(aliases),'duplicate paths in',dest)
if __name__=='__main__':main()
