from pathlib import Path
import csv,json,gzip,openpyxl,re
P=Path('/data/CARBON_PROJECT');R=P/'revision_v20'
B=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/evidence_expansion_r4/external_asb')
with (P/'revision_v4/work_variants/EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv').open(encoding='utf8') as f:d=list(csv.DictReader(f))
ev=json.loads(gzip.decompress((P/'revision_v19/work_fig5/frozen_inputs/work_variants/source_inputs/EVENTS_INCLUDED.json.gz').read_bytes())); em={e['event_id']:e for e in ev}
# Published Figure4e exact pooled genomic counts, complemented only to source RNA orientation.
targets={'MAP2K3':{'BUD13':('C','G',24,0),'PRPF4':('C','G',26,0),'SF3A3':('C','G',30,0),'SF3B4':('C','G',26,0)},'PDE4DIP':{'BUD13':('C','U',42,15),'MTPAP':('C','U',24,5),'ZNF622':('C','U',46,12)},'LARS':{'U2AF2':('C','U',86,39)},'OTUB1':{'BUD13':('C','G',28,10),'PRPF8':('C','G',36,14),'SF3B4':('C','G',21,4)},'IL17RB':{'SF3B4':('C','U',46,15)}}
comp={'A':'T','T':'A','C':'G','G':'C'}
matched=[]
for gene,t in targets.items():
    for rb,(a,b,na,nb) in t.items():
        cand=[r for r in d if gene in r['gene_names'].split(';') and r['RBP']==rb]
        hits=[]
        for r in cand:
            e=em[r['event_id']]; raw=r['source_raw_counts'];counts={base:sum(int(n) for bb,n in re.findall(r'([ACGT]):(\d+)',raw) if bb==base) for base in 'ACGT'}
            # Source strand refers to hg19 + genomic base counts.
            ra=a.replace('U','T');rbv=b.replace('U','T')
            if e['source_strand']=='-':ra,rbv=comp[ra],comp[rbv]
            if counts[ra]==na and counts[rbv]==nb:hits.append(r)
        print(gene,rb,'exact count matched',[(r['event_id'],r['variant_key'],r['strand'],r['A_native_credit'],r['S_credit'],r['A_flat_credit'],r['A_permmean_credit']) for r in hits])
        for r in hits:matched.append(dict(r,published_gene=gene,published_RBP=rb,published_preferred_rna=a,published_other_rna=b,source_count_match='exact pooled counts from Figure4e; primer verification required'))
with (R/'analysis/FIG4E_COUNT_MATCHED_EVENTS.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(matched[0]));w.writeheader();w.writerows(matched)
wb=openpyxl.load_workbook(B/'41467_2019_9292_MOESM4_ESM.xlsx',read_only=True,data_only=True)
print('PRIMERS',[(i,row) for i,row in enumerate(wb.active.iter_rows(values_only=True),1) if any(g in str(row) for g in ['LARS','MAP2K3','PDE4DIP','OTUB1','IL17RB','DNAJB2','CERS2','UBE2I','PDXK','ACVR1B','SLC39A13'])])
# keep terse details of all locus events for exact main candidate.
for key in ['chr1:148974518:C:T','chr11:63996830:G:C']:
    print(key,[(r['RBP'],r['cell'],r['ref_preferred'],r['A_native_credit'],r['S_credit']) for r in d if r['variant_key']==key])
