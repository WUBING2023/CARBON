from pathlib import Path
import csv,json,gzip,openpyxl
P=Path('/data/CARBON_PROJECT');R=P/'revision_v20'
B=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/evidence_expansion_r4/external_asb')
with (P/'revision_v4/work_variants/EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv').open(encoding='utf8') as f:rows=list(csv.DictReader(f))
events=json.loads(gzip.decompress((P/'revision_v19/work_fig5/frozen_inputs/work_variants/source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()));em={e['event_id']:e for e in events}
wb=openpyxl.load_workbook(B/'41467_2019_9292_MOESM4_ESM.xlsx',read_only=True,data_only=True)
prim=list(wb.active.iter_rows(values_only=True))
for i in range(21,28):print('primer',i,prim[i-1])
pairs={'MAP2K3':(6,8),'PDE4DIP':(12,14),'LARS':(18,20),'OTUB1':(24,26),'IL17RB':(30,32),'DNAJB2':(36,38),'CERS2':(42,44),'SLC39A13':(48,50),'ACVR1B':(57,60),'PDXK':(63,66),'UBE2I':(69,69)}
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
out=[];fail=[]
for gene,indices in pairs.items():
    seqs=[''.join(str(prim[i-1][1]).upper().split()) for i in indices]
    if gene=='UBE2I':
        # Published Figure5e is C/G. Only the G mutagenesis primer is supplied;
        # Change the focal G at index11 to C, verified against the Figure5e
        # source hg19 coordinate chr16:1374951 and both frozen alleles.
        assert seqs[0]=='GTGTGCGGTCTGGATTCGC'
        seqs[1]=seqs[0][:11]+'C'+seqs[0][12:]
    if len(seqs[0])!=len(seqs[1]):
        fail.append({'gene':gene,'reason':'Unequal-length primers; inspect alternative allele or Figure5 coordinate','primers':seqs});continue
    diff=[i for i,(a,b) in enumerate(zip(*seqs)) if a!=b]
    if len(diff)!=1:
        fail.append({'gene':gene,'reason':'not single-difference primer pair','primers':seqs});continue
    hits=[]
    for r in rows:
        if gene not in r['gene_names'].split(';'):continue
        c=em[r['event_id']]['context'];ref=c['reference_sequence_rna_direction_dna'];other=c['minor_sequence_dna'] if c['reference_rna_base_dna']==c['major_rna_base_dna'] else c['major_sequence_dna']
        for orientation,ss,ix in [('direct',seqs,diff[0]),('reverse_complement',[rc(s) for s in seqs],len(seqs[0])-diff[0]-1)]:
            st=50-ix
            if st<0 or st+len(ss[0])>101:continue
            if {ref[st:st+len(ss[0])],other[st:st+len(ss[0])]}==set(ss):
                hits.append(dict(r,published_gene=gene,primer_rows=';'.join(map(str,indices)),primer_orientation_in_model_input=orientation,primer_reference_index=ss.index(ref[st:st+len(ss[0])]),primer_snv_offset=ix,primer_start_in101=st,isolated_reporter_snv_exact_match=True,primer_inferred_alternative=gene=='UBE2I'))
    keys=sorted({r['variant_key'] for r in hits})
    print(gene,'EXACT',keys,[(r['RBP'],r['A_native_credit'],r['S_credit'],r['A_flat_credit'],r['A_permmean_credit']) for r in hits])
    if len(keys)>1:raise ValueError((gene,keys))
    out.extend(hits)
    if not hits:fail.append({'gene':gene,'reason':'no exact centered primer match in frozen qualified inputs','primers':seqs})
with (R/'analysis/PUBLISHED_REPORTER_EXACT_MATCHES.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
(R/'analysis/REPORTER_MATCH_FAILURES.json').write_text(json.dumps(fail,indent=2),encoding='utf8')
print('FAILED',fail)
