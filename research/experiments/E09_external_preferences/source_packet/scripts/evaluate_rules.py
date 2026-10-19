from literal_tables import read_csv_literal
"""Evaluate locked rule extraction against genuine external RBNS/RNAcompete tables."""
import json, pathlib, itertools, hashlib, tarfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=pathlib.Path(__file__).resolve().parents[1]
KMERS=[''.join(k) for k in itertools.product('ACGT',repeat=5)]
ARMS=['aligned','sequence_only','mean','permuted']

def ranking(score,eligible=None):
    keys=KMERS if eligible is None else eligible
    return sorted(keys,key=lambda k:(-float(score[k]),k))

def truth():
    res={};metadata=[]
    for protein,accession in [('HNRNPK','ENCFF715DIL'),('RBFOX2','ENCFF002DFE')]:
        df=read_csv_literal(ROOT/'data'/f'{accession}.tsv',sep='\t').rename(columns=lambda s:str(s).strip())
        keys=df.iloc[:,0].astype(str).str.upper().str.replace('U','T')
        assert len(keys)==1024 and set(keys)==set(KMERS) and not keys.duplicated().any()
        vals=df.iloc[:,1:].astype(float);assert np.isfinite(vals.to_numpy()).all()
        chosen=vals.max().idxmax()
        for c in vals:
            name=f'{accession}:{c}'
            res.setdefault(protein,{})[name]={k:float(v) for k,v in zip(keys,vals[c])}
            metadata.append(dict(protein=protein,assay='RBNS',experiment=accession,column=c,label=name,primary=(c==chosen),n_kmers=1024,min=float(vals[c].min()),max=float(vals[c].max()),U_to_T=True))
    file=ROOT/'sources'/'Zscores_420_origrncmpt.txt'
    with file.open() as f:header=f.readline().split()[1:]
    cols=[0]+[header.index(v)+1 for v in ['RNCMPT00268','RNCMPT00269']]
    df=read_csv_literal(file,sep=r'\s+',skiprows=1,header=None,usecols=cols)
    keys=df.iloc[:,0].astype(str).str.upper().str.replace('U','T')
    expected=set(''.join(k) for k in itertools.product('ACGT',repeat=7))
    assert len(keys)==16382 and not keys.duplicated().any() and set(keys).issubset(expected)
    missing=sorted(expected-set(keys));assert missing==['GAAGAGC','GCTCTTC']
    for jj,experiment in enumerate(['RNCMPT00268','RNCMPT00269'],start=1):
        sums={k:0. for k in KMERS};nums={k:0 for k in KMERS}
        vals=df.iloc[:,jj].to_numpy(float);assert np.isfinite(vals).all()
        for k,z in zip(keys,vals):
            for offset in range(3):
                kk=k[offset:offset+5];sums[kk]+=float(z);nums[kk]+=1
        # Repeated 5-mers within one 7-mer count repeatedly, matching RBPNet.
        assert min(nums.values())>0
        res.setdefault('PTBP1',{})[experiment]={k:sums[k]/nums[k] for k in KMERS}
        metadata.append(dict(protein='PTBP1',assay='RNAcompete',experiment=experiment,column='7merZ_occurrence_mean',label=experiment,primary=True,n_7mers=len(keys),n_kmers=1024,occurrences_per_5mer_min=min(nums.values()),occurrences_per_5mer_max=max(nums.values()),missing_7mers=';'.join(missing),missing_imputation=False,U_to_T=True))
    pd.DataFrame(metadata).to_csv(ROOT/'results'/'EXTERNAL_COLUMNS.csv',index=False)
    rows=[dict(protein=p,label=l,kmer=k,score=v) for p,labels in res.items() for l,scores in labels.items() for k,v in scores.items()]
    pd.DataFrame(rows).to_csv(ROOT/'results'/'EXTERNAL_5MER_SCORES.csv',index=False)
    return res,{m['label'] for m in metadata if m['primary']}

def evaluate():
    extract=ROOT/'results'/'extracted';extract.mkdir(exist_ok=True)
    if not (extract/'outputs').exists():
        with tarfile.open(ROOT/'results'/'ig_outputs.tar.gz') as t:
            for member in t:
                target=(extract/member.name).resolve();assert target.is_relative_to(extract.resolve())
                t.extract(member,extract)
    truths,primary=truth();metrics=[];kmerranks=[];audits=[]
    output=extract/'outputs'
    for task in sorted(output.iterdir()):
        protein,cell=task.name.split('_')[:2];native=None;sequences=None
        for seed in [11,23,47]:
            for arm in ARMS:
                p=task/f'seed{seed}'/arm
                j=json.loads((p/'rules.json').read_text());m=json.loads((p/'summary.json').read_text())
                a=np.load(p/'attributions.npz');seq=a['sequence'];keys=j['kmers']
                assert m['n_selected']==1024 and len(seq)==1024 and set(keys)==set(KMERS)
                if native is None:native=seq;sequences=a
                assert np.array_equal(seq,native)
                source_chrom=a['chromosome'];source_row=a['source_row'];phi=a['conservation']
                transformed_phi=phi.copy()
                if arm=='mean':transformed_phi[:]=phi.mean(1,keepdims=True)
                elif arm=='permuted':
                    for ii,(c,r) in enumerate(zip(source_chrom,source_row)):
                        b=hashlib.sha256(f'20260929|{task.name}|{int(c)}|{int(r)}'.encode()).digest()
                        perm=np.random.default_rng(int.from_bytes(b[:8],'big')).permutation(101)
                        transformed_phi[ii]=phi[ii,perm]
                audits.append(dict(task=task.name,seed=seed,arm=arm,n_selected=m['n_selected'],n_total=m['n_total'],steps=m['steps'],scalar=m['scalar'],checkpoint_sha256=m['checkpoint_sha256'],input_sha256=m['input_sha256'],sequence_sha256=m['selected_sequence_sha256'],native_phi_sha256=hashlib.sha256(phi.tobytes()).hexdigest(),conditional_phi_sha256=hashlib.sha256(transformed_phi.tobytes()).hexdigest() if arm!='sequence_only' else None,ig_phi_baseline_equals_input=True,completeness_error=m['mean_abs_completeness_error'],mean_abs_output_difference=m['mean_abs_output_difference'],elapsed_seconds=m['elapsed_seconds']))
                mean=dict(zip(keys,j['mean_occurrence']));winner=dict(zip(keys,j['winner_sum']));selected=[k for k,n in zip(keys,j['winner_n']) if n>0];top=ranking(winner,selected)[:20]
                for label,actual in truths[protein].items():
                    gt=ranking(actual)[:20];recall=len(set(top)&set(gt))/20;rho=float(spearmanr([mean[k] for k in KMERS],[actual[k] for k in KMERS]).statistic)
                    metrics.append(dict(task=task.name,protein=protein,cell=cell,seed=seed,arm=arm,label=label,primary_column=label in primary,recall20=recall,spearman_all5mer=rho,n_windows=1024,n_kmers=1024,n_recovered_primary=len(selected),matched_training=True))
                wrank={k:i+1 for i,k in enumerate(ranking(winner,selected))}
                for k in KMERS:kmerranks.append(dict(task=task.name,protein=protein,cell=cell,seed=seed,arm=arm,kmer=k,mean_occurrence_ig=mean[k],winner_sum_ig=winner[k],winner_n=j['winner_n'][keys.index(k)],primary_rank=wrank.get(k,np.nan)))
        # Baselines operate on exactly the native, arm-common sequences.
        seq=native;phi=sequences['conservation'];freq={k:0. for k in KMERS};pm={k:0. for k in KMERS};pn={k:0 for k in KMERS};pw={k:0. for k in KMERS};pwn={k:0 for k in KMERS}
        for s,ph in zip(seq,phi):
            v=np.lib.stride_tricks.sliding_window_view(ph,5).mean(1);best=int(np.argmax(v));pw[s[best:best+5]]+=float(v[best]);pwn[s[best:best+5]]+=1
            for i in range(97):k=s[i:i+5];freq[k]+=1;pm[k]+=float(v[i]);pn[k]+=1
        pm={k:pm[k]/pn[k] for k in KMERS};consensus={k:consensus_score(protein,k) for k in KMERS}
        for arm,sec,pri,eligible in [('frequency',freq,freq,KMERS),('phi_priority',pm,pw,[k for k in KMERS if pwn[k]>0]),('known_consensus',consensus,consensus,KMERS)]:
            top=ranking(pri,eligible)[:20]
            for label,actual in truths[protein].items():
                metrics.append(dict(task=task.name,protein=protein,cell=cell,seed=0,arm=arm,label=label,primary_column=label in primary,recall20=len(set(top)&set(ranking(actual)[:20]))/20,spearman_all5mer=float(spearmanr([sec[k] for k in KMERS],[actual[k] for k in KMERS]).statistic),n_windows=1024,n_kmers=1024,n_recovered_primary=len(eligible),matched_training=False))
    pd.DataFrame(audits).to_csv(ROOT/'results'/'IG_AUDIT.csv',index=False)
    pd.DataFrame(kmerranks).to_csv(ROOT/'results'/'MODEL_5MER_SCORES.csv',index=False)
    df=pd.DataFrame(metrics);df.to_csv(ROOT/'results'/'ALL_METRICS.csv',index=False)
    pri=df[df.primary_column]
    # RNAcompete experiments are averaged, never selected for model agreement.
    bytask=pri.groupby(['protein','cell','arm'],as_index=False)[['recall20','spearman_all5mer']].mean();bytask.to_csv(ROOT/'results'/'BACKGROUND_METRICS.csv',index=False)
    byprotein=bytask.groupby(['protein','arm'],as_index=False)[['recall20','spearman_all5mer']].mean();byprotein.to_csv(ROOT/'results'/'PROTEIN_METRICS.csv',index=False)
    print(byprotein.to_string(index=False))
    comparisons=[]
    for protein in byprotein.protein.unique():
        sub=byprotein[byprotein.protein==protein].set_index('arm')
        for control in ['sequence_only','mean','permuted','frequency','phi_priority','known_consensus']:
            comparisons.append(dict(protein=protein,control=control,delta_recall20=float(sub.loc['aligned','recall20']-sub.loc[control,'recall20']),delta_spearman=float(sub.loc['aligned','spearman_all5mer']-sub.loc[control,'spearman_all5mer'])))
    pd.DataFrame(comparisons).to_csv(ROOT/'results'/'PAIRED_PROTEIN_DIFFERENCES.csv',index=False)

def consensus_score(protein,k):
    if protein=='HNRNPK':return k.count('C')
    if protein=='PTBP1':return k.count('C')+k.count('T')
    return max([len(k[i:j]) for i in range(5) for j in range(i+1,6) if k[i:j] in 'TGCATG']+[0])

if __name__=='__main__':evaluate()
