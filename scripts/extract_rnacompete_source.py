#!/usr/bin/env python3
"""Verify the pinned original measured RNAcompete matrix, or extract it outside Git."""
from pathlib import Path
import argparse,hashlib,json,tarfile
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'paper/source_data/S23_external_order_preferences/original/sources/Zscores_420_origrncmpt.tar.gz'
ARCHIVE_SHA='bf65e4bdd12fb1f8fd5db585e5165db24cc1dba6a1b2708f2cc0e12255a676bd'
MATRIX_SHA='c4225eefb9561dfdb360858aad1a50545176df23638e62b077b7eb683877d6bd'
MATRIX_NAME='Zscores_420_origrncmpt.txt'
def checksum(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,default=ARCHIVE);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--verify-only',action='store_true');g.add_argument('--destination',type=Path);a=p.parse_args()
 assert a.archive.stat().st_size==51769595,'Unexpected original archive size'
 assert checksum(a.archive)==ARCHIVE_SHA,'Original archive SHA256 mismatch'
 target=None
 if a.destination:
  target=a.destination.resolve()
  if target==ROOT or ROOT in target.parents:raise SystemExit('Extract the 115 MB matrix outside the Git repository.')
  if target.exists():raise SystemExit('Destination exists; choose a new file.')
  target.parent.mkdir(parents=True,exist_ok=True)
 h=hashlib.sha256();size=0;lines=0;output=target.open('xb') if target else None
 try:
  with tarfile.open(a.archive,'r:gz') as t:
   member=t.getmember(MATRIX_NAME);assert member.isfile() and member.size==115308739
   with t.extractfile(member) as f:
    for b in iter(lambda:f.read(1024*1024),b''):
     h.update(b);size+=len(b);lines+=b.count(b'\n')
     if output:output.write(b)
 finally:
  if output:output.close()
 assert size==115308739 and h.hexdigest()==MATRIX_SHA,'Extracted measured matrix identity mismatch'
 print(json.dumps(dict(passed=True,archive_sha256=ARCHIVE_SHA,matrix_sha256=h.hexdigest(),matrix_bytes=size,newline_count=lines,extracted=bool(target),scope='Original measured source byte identity only; no 7-mer to 5-mer mapping or model rerun.'),indent=2))
if __name__=='__main__':main()
