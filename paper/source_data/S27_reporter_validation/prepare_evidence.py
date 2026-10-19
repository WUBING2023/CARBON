from pathlib import Path
import csv,json,gzip,hashlib,shutil,openpyxl
P=Path('/data/CARBON_PROJECT');R=P/'revision_v20';V=P/'revision_v19'
B=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/evidence_expansion_r4/external_asb')
for d in ['manuscript','work_fig5','work_disease','build']:
    if not (R/d).exists():shutil.copytree(V/d,R/d,ignore=shutil.ignore_patterns('__pycache__','pdf_qa','*.log','*.aux','*.fls','*.fdb_latexmk','*.out','*.synctex.gz'))
(R/'qa').mkdir(exist_ok=True)
F=R/'work_fig5/frozen_inputs/reporter_validation';F.mkdir(exist_ok=True)
for name in ['41467_2019_9292_MOESM4_ESM.xlsx','41467_2019_9292_MOESM6_ESM.xlsx']:
    shutil.copyfile(B/name,F/name)
shutil.copyfile(R/'analysis/PUBLISHED_REPORTER_EXACT_MATCHES.csv',F/'ALL24_REPORTER_LOCUS_ASB_EVENTS.csv')
def read(p):
    with p.open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))
ev=read(F/'ALL24_REPORTER_LOCUS_ASB_EVENTS.csv');assert len(ev)==24 and len({r['variant_key'] for r in ev})==11
wb=openpyxl.load_workbook(F/'41467_2019_9292_MOESM6_ESM.xlsx',read_only=True,data_only=True)
meas=[];sites=[]
for sheet in ['Figure 4e','Figure 5e','Figure 6c']:
    for i,row in enumerate(wb[sheet].iter_rows(values_only=True),1):
        gene=str(row[0]);gene='SLC39A13' if gene=='SCL39A13' else gene
        if gene not in {r['published_gene'] for r in ev}:continue
        splice=sheet=='Figure 4e' or (sheet=='Figure 6c' and gene in ['DNAJB2','CERS2'])
        values=row[2:8] if splice else row[1:7];bases=row[8:14] if splice else row[7:13]
        endpoint=('exon '+str(row[1])+' (%)' if gene in ['MAP2K3','PDE4DIP','LARS','CERS2'] else 'intron '+str(row[1])+' (%)') if splice else 'mCherry/eYFP normalized expression'
        matches=[r for r in ev if r['published_gene']==gene];keys={r['variant_key'] for r in matches};assert len(keys)==1
        for j,(val,base) in enumerate(zip(values,bases)):
            meas.append(dict(gene=gene,variant_key=next(iter(keys)),source_sheet=sheet,source_excel_row=i,source_value_column=j+(3 if splice else 2),RNA_allele=str(base).strip("'"),biological_replicate=j%3+1,endpoint=endpoint,value=float(val),assay_cell='HeLa',source_units='percent' if splice else 'normalized_ratio'))
        sites.append(dict(gene=gene,variant_key=next(iter(keys)),assay=endpoint,source_sheet=sheet,source_excel_row=i,assay_cell='HeLa',biological_replicates_per_allele=3,ASB_events=len(matches),split=matches[0]['split'],Full_correct=sum(float(r['A_native_credit']) for r in matches),Sequence_correct=sum(float(r['S_credit']) for r in matches)))
assert len(meas)==66 and len(sites)==11
def write(p,rows):
    with p.open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write(F/'ALL66_REPORTER_BIOLOGICAL_MEASUREMENTS.csv',meas);write(F/'ALL11_REPORTER_VARIANTS.csv',sites)
events=json.loads(gzip.decompress((R/'work_fig5/frozen_inputs/work_variants/source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()));em={e['event_id']:e for e in events}
selected=['K562:UPF1:chr12_52388891_+','K562:GEMIN5:chr1_150940625_-']
packet={eid:{'event':em[eid],'prediction':next(r for r in ev if r['event_id']==eid),'all_locus_events':[r for r in ev if r['variant_key']==next(q['variant_key'] for q in ev if q['event_id']==eid)],'reporter':[r for r in meas if r['gene']==next(q['published_gene'] for q in ev if q['event_id']==eid)]} for eid in selected}
(F/'SELECTED_CASES.json').write_text(json.dumps(packet,indent=2),encoding='utf8')
protocol=json.loads((R/'analysis/CASE_SELECTION_PROTOCOL.json').read_text('utf8'))
protocol['reporter_addendum']={'source':'Yang et al.2019 Nature Communications10:1338, doi:10.1038/s41467-019-09292-w','availability':'All11 two-allele reporter loci in original Figures4e/5e/6c; all24 qualified existing ASB events at exact loci;66 biological measurements. No new predictions or assays.','matching':'Both mutagenesis primers exactly match reference/alternate101 windows including focal position50, with direct/reverse-complement orientation recorded. UBE2I alternative C primer inferred from published C/G and verified against original hg19chr16:1374951 and both frozen allele windows; only G primer supplied. Source typo SCL39A13 is normalized explicitly to SLC39A13.','main_selection':'ACVR1B selected posthoc for a single-RBP expression-reporter locus where Full agrees with ASB and Sequence/flat/permutation mean disagree; individual perm2 also agrees and is displayed. CERS2 selected for exact GWAS rs267738 plus published isolated-allele splicing and GEMIN5-KD support; Full and Sequence agree for GEMIN5, and complete locus YBX3 Full failure displayed. No selected-case success estimate.','not_established':'CARBON does not predict the reporter expression/splicing measurement and no direct quantitative model-to-phenotype endpoint or RBP-mediated causality test is computed. All source ASB/model labels exposed; training chromosomes remain labelled.'}
(F/'PROTOCOL_AND_SELECTION.json').write_text(json.dumps(protocol,indent=2),encoding='utf8')
for n in ['CASE_CENSUS.json','ALL_VARIANT_CELL_STRAND_GROUPS.csv','FULL_ONLY_STABLE_EVENT_ILLUSTRATIONS.csv']:
    shutil.copyfile(R/'analysis'/n,F/n)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(F/'SOURCE_HASHES.json').write_text(json.dumps([{'file':f.name,'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(F.iterdir()) if f.is_file()],indent=2),encoding='utf8')
print(json.dumps({'reporter_loci':11,'all_ASB_events':24,'assay_values':66,'Full_correct':sum(float(r['A_native_credit']) for r in ev),'Sequence_correct':sum(float(r['S_credit']) for r in ev),'sites':sites},indent=2))
