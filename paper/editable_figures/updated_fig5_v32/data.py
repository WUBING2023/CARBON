"""Native-vector Figure5 redesign from immutable, frozen inputs only."""
from pathlib import Path
import csv, json, gzip, hashlib, shutil, os
os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, AsinhNorm
import fitz

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PROJECT=ROOT.parent
V4=PROJECT/'revision_v4'
R=HERE/'frozen_inputs/work_variants'
S=HERE/'frozen_inputs/S23_external_order_preferences'
P=np.load(R/'b1_reference/INTERVENTION_PACKET.npz')
W,H=178.,162.
FULL,SEQ,FREQ,FLAT,PERM='#E64B35','#4DBBD5','#3C5488','#00A087','#8491B4'
INK,GRAY,PALE='#263238','#74808B','#E1E6EA'
BASE=json.loads((HERE/'frozen_inputs/PALETTE_V3.json').read_text('utf8'))['nucleotides']
ARMS=['A_native','S','A_flat','A_permmean','ref']
LABELS=['Full','Sequence','Flat Φ','Perm. Φ','REF rule']
COLORS=[FULL,SEQ,FLAT,PERM,FREQ]
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,
 'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.linewidth':.55,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})

def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def writecsv(p,data):
    with p.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

metrics=rows(S/'FROZEN_PRIMARY_METRICS.csv')
contrasts=rows(S/'FROZEN_PRIMARY_CONTRASTS.csv')
proteins=P['cohort_proteins'].tolist()
lookup={(x['protein'],x['arm']):x for x in metrics}
assert len(proteins)==25 and len(metrics)==75
example=rows(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv')
events=json.loads(gzip.decompress((R/'source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()))
ev=next(e for e in events if e['event_id']==example[0]['event_id'])
ref=ev['context']['reference_sequence_rna_direction_dna'].replace('T','U')
alt=ev['context']['minor_sequence_dna'].replace('T','U')
phi=np.asarray(ev['context']['phi_raw'])
assert len(example)==10 and ref[50]=='C' and alt[50]=='U'
cohort=rows(R/'STRICT396_WITH_ANNOTATIONS.csv')
cohort.sort(key=lambda x:(float(x['phi_center']),int(x['chrom'][3:]),int(x['pos1']),x['ref'],x['alt'],x['cell'],x['RBP'],x['event_id']))
assert len(cohort)==396 and len({x['variant_key'] for x in cohort})==287
assert len({x['RBP'] for x in cohort})==85 and len({x['cluster_id'] for x in cohort})==119
quart=[r for r in rows(R/'DIRECTION_BY_CONSERVATION.csv') if r['cohort']=='strict_test']
summary={r['arm']:r for r in rows(R/'DIRECTION_SUMMARY.csv') if r['cohort']=='strict_test'}
two={r['arm']:r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv') if r['cohort']=='strict_test'}
assoc=[r for r in rows(R/'CONSERVATION_ASSOCIATIONS.csv') if r['cohort']=='strict_test' and r['arm']=='A_native']
assign={r['event_id']:int(r['phi_quartile']) for r in rows(R/'QUARTILE_ASSIGNMENTS_strict_test.csv')}
counts=[sum(assign[r['event_id']]==i for r in cohort) for i in range(1,5)]
assert counts==[116,101,86,93]
