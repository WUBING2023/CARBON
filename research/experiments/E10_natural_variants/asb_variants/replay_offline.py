"""Recompute archived endpoints/annotations/figures locally, without downloads."""
from pathlib import Path
import subprocess,sys,argparse,json
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--metrics-only',action='store_true');args=p.parse_args()
assert (R/'source_inputs/EVENTS_INCLUDED.json.gz').exists(), 'Copy source_inputs with this package first.'
assert (R/'annotation_ncbi/LOOKUP_ALL_VARIANTS.json').exists(), 'Copy annotation caches first.'
scripts=['analyze_variants.py','verify_weighted_refits.py','parse_annotations.py','parse_clinvar_cache.py','finalize_annotations.py','verify_region_sensitivity.py','audit_sources.py']
if not args.metrics_only:scripts+=['render_variants.py','render_supplement.py','export_continued.py']
for s in scripts:
    print('RUN '+s,flush=True);subprocess.run([sys.executable,str(R/s)],cwd=R,check=True)
print('Offline replay completed; original frozen protocol and acquisition logs retained.')
