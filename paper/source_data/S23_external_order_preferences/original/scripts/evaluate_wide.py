from literal_tables import read_csv_literal
"""Availability-frozen28-protein evaluation; no outcome-based selection."""
import pathlib,json,itertools,tarfile,hashlib
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from evaluate_rules import ROOT,KMERS
from evaluate_controlled import pairs,concordance

ARMS=['full','sequence','frozen_flat','frozen_perm']
def external(manifest):
    mapping={}
    for spec in manifest:
        for experiment in spec['experiments']:mapping[experiment]=spec['protein']
    p=ROOT/'sources'/'Zscores_420_origrncmpt.txt'
    with p.open() as f:header=f.readline().split()[1:]
    labs=sorted(mapping);assert all(k in header for k in labs)
    df=read_csv_literal(p,sep=r'\s+',skiprows=1,header=None,usecols=[0]+[header.index(k)+1 for k in labs]);columns=[header[int(c)-1] for c in df.columns if int(c)>0];words=df.iloc[:,0].str.replace('U','T');assert len(words)==16382 and len(set(words))==16382
    true={};rows=[];meta=[]
    author=read_csv_literal(ROOT/'sources'/'S1.In_Progress.v6.For_all.DR-2.csv',sep='\t',keep_default_na=False)
    qa=read_csv_literal(ROOT/'sources'/'Fullset_z_scores-setA-vs-setB-comparison.txt',sep=r'\s+',skiprows=1,header=None,names=['experiment','setA_B_pearson','top10_overlap','top100_overlap','top1000_overlap']).set_index('experiment')
    for jj,label in enumerate(columns,1):
        val=df.iloc[:,jj].to_numpy(float);assert np.isfinite(val).all();summ=np.zeros(1024);n=np.zeros(1024,dtype=int);ki={k:i for i,k in enumerate(KMERS)}
        for word,z in zip(words,val):
            for offset in range(3):i=ki[word[offset:offset+5]];summ[i]+=float(z);n[i]+=1
        assert (n>0).all();scores=summ/n;true.setdefault(mapping[label],{})[label]=scores
        row=author[author.Motif_ID.str.split(',').apply(lambda x:label in x)].iloc[0]
        meta.append(dict(protein=mapping[label],experiment=label,human_gene_id=row.Gene_systematic_ID,source=row.Source,source_reference=row.Source_Reference,project=row.Project,domains=row.Domains,construct_aa_length=len(row.Construct_AA_seq),hyb_id=row.Hyb_ID,setA_B_pearson=float(qa.loc[label,'setA_B_pearson']) if label in qa.index else np.nan,n_7mers=16382,n_5mers=1024,min_5mer_occurrences=int(n.min()),max_5mer_occurrences=int(n.max()),all_measured=True,exposed_protein=mapping[label] in ['HNRNPK','PTBP1','RBFOX2']))
        for k,v in zip(KMERS,scores):rows.append(dict(protein=mapping[label],experiment=label,kmer=k,z5=float(v),n_occurrences=int(n[ki[k]])))
    pd.DataFrame(rows).to_csv(ROOT/'results'/'WIDE_EXTERNAL_5MER_SCORES.csv',index=False);pd.DataFrame(meta).to_csv(ROOT/'results'/'WIDE_EXPERIMENT_METADATA.csv',index=False)
    return true

def run():
    man=json.loads((ROOT/'WIDE_MANIFEST.json').read_text());true=external(man);folder=ROOT/'results'/'wide_extracted';folder.mkdir(exist_ok=True)
    with tarfile.open(ROOT/'results'/'wide_outputs.tar.gz') as tf:
        for m in tf:
            assert (folder/m.name).resolve().is_relative_to(folder.resolve());tf.extract(m,folder)
    pp,groups=pairs();metrics=[];scores=[];audit=[];domain=[];restricted=[]
    for spec in man:
        task=spec['task'];protein=spec['protein'];cell=task.split('_')[1];x=np.load(folder/'wide_inputs'/f'{task}.npz');words=x['sequence'];phi=x['x'][:,:,4];freq=np.zeros(1024);phiv=np.zeros(1024);nn=np.zeros(1024);ki={k:i for i,k in enumerate(KMERS)}
        for ss,p in zip(words,phi):
            vals=np.lib.stride_tricks.sliding_window_view(p,5).mean(1)
            for off in range(97):i=ki[ss[off:off+5]];freq[i]+=1;phiv[i]+=float(vals[off]);nn[i]+=1
        valid=nn>0;phiv[valid]/=nn[valid];phiv[~valid]=np.nan;preds={'frequency':freq,'phi_only':np.zeros(1024)}
        domain.append(dict(task=task,protein=protein,n_native_source_windows=len(words),n_observed_words=int(valid.sum()),n_missing_words=int((~valid).sum()),missing_words=';'.join(k for k,g in zip(KMERS,valid) if not g),missing_phi_not_imputed=True))
        for base in 'ACGT':preds[base+'_count']=np.array([k.count(base) for k in KMERS],float)
        for arm in ARMS:
            path=folder/'wide_outputs'/task/arm;arr=np.load(path/'responses.npz');meta=json.loads((path/'summary.json').read_text());assert list(arr['kmers'])==KMERS and arr['delta_log_count'].shape==(16,1024);assert np.array_equal(arr['sequence'],words[:16]);assert meta['checkpoint_sha256']==spec['sequence_checkpoint_sha256' if arm=='sequence' else 'full_checkpoint_sha256'];preds[arm]=arr['delta_log_count'].mean(0)
            cp=arr['conditional_conservation'];orig=arr['native_conservation']
            if arm=='full':assert np.array_equal(cp,orig)
            elif arm=='sequence':assert cp.shape==(0,101)
            elif arm=='frozen_flat':assert np.allclose(cp,orig.mean(1,keepdims=True))
            elif arm=='frozen_perm':assert np.array_equal(np.sort(cp,1),np.sort(orig,1))
            audit.append(dict(task=task,protein=protein,arm=arm,checkpoint_hash_pass=True,input_identity_pass=True,phi_transform_pass=True,n_contexts=16,n_kmers=1024,elapsed_seconds=meta['elapsed_seconds'],exposed_protein=spec['exposed_protein']))
        for arm,pred in preds.items():
            for k,v in zip(KMERS,pred):scores.append(dict(task=task,protein=protein,cell=cell,arm=arm,kmer=k,score=float(v),exposed_protein=spec['exposed_protein']))
            for lab,actual in true[protein].items():
                cc,n,t,e=concordance(pred,actual,pp);rho=float(spearmanr(pred,actual).statistic) if np.ptp(pred)>0 else np.nan
                metrics.append(dict(task=task,protein=protein,cell=cell,arm=arm,experiment=lab,spearman_all5mer=rho,within_composition_concordance=cc,n_pairs=n,predicted_ties=t,experimental_ties=e,n_contexts=16,n_kmers=1024,exposed_protein=spec['exposed_protein'],matched_training_causal_comparison=False))
        observed_pairs=pp[valid[pp[:,0]]&valid[pp[:,1]]]
        for arm,pred in dict(preds,phi_association=phiv).items():
            for lab,actual in true[protein].items():
                cc,n,t,e=concordance(pred,actual,observed_pairs);rho=float(spearmanr(pred[valid],actual[valid]).statistic) if np.ptp(pred[valid])>0 else np.nan
                restricted.append(dict(task=task,protein=protein,cell=cell,arm=arm,experiment=lab,spearman_observed_word_domain=rho,within_composition_concordance=cc,n_words=int(valid.sum()),n_pairs=n,scope='same observed native5mer domain for every method',exposed_protein=spec['exposed_protein']))
    assert len(audit)==168 and len({a['task'] for a in audit})==42 and len({a['protein'] for a in audit})==28
    df=pd.DataFrame(metrics);df.to_csv(ROOT/'results'/'WIDE_ALL_METRICS.csv',index=False);pd.DataFrame(scores).to_csv(ROOT/'results'/'WIDE_ALL_5MER_PREDICTIONS.csv',index=False);pd.DataFrame(audit).to_csv(ROOT/'results'/'WIDE_EXECUTION_QA.csv',index=False);pd.DataFrame(domain).to_csv(ROOT/'results'/'WIDE_OBSERVED_DOMAIN_QA.csv',index=False);pd.DataFrame(restricted).to_csv(ROOT/'results'/'WIDE_NATIVE_PHI_COMMON_DOMAIN.csv',index=False)
    bg=df.groupby(['protein','cell','arm','exposed_protein'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();bg.to_csv(ROOT/'results'/'WIDE_BACKGROUND_METRICS.csv',index=False)
    prot=bg.groupby(['protein','arm','exposed_protein'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();prot.to_csv(ROOT/'results'/'WIDE_PROTEIN_METRICS.csv',index=False)
    rng=np.random.default_rng(20261004);summary=[];diffrows=[]
    for name,sub in [('all28',prot),('new25',prot[~prot.exposed_protein]),('exposed3',prot[prot.exposed_protein])]:
        proteins=sorted(sub.protein.unique());n=len(proteins);ix=rng.integers(0,n,(10000,n))
        for metric in ['spearman_all5mer','within_composition_concordance']:
            tab=sub.pivot(index='protein',columns='arm',values=metric).loc[proteins]
            for arm in tab:
                vals=tab[arm].to_numpy();b=vals[ix].mean(1)
                summary.append(dict(subset=name,endpoint=metric,arm=arm,n_proteins=n,mean=float(vals.mean()),median=float(np.median(vals)),protein_bootstrap_low=float(np.quantile(b,.025)),protein_bootstrap_high=float(np.quantile(b,.975)),n_proteins_above_zero_or_half=int((vals>(0 if metric=='spearman_all5mer' else .5)).sum())))
            for arm in ['sequence','frozen_flat','frozen_perm','frequency','phi_only','A_count','C_count','G_count','T_count']:
                if not np.isfinite(tab[arm]).all():continue
                d=(tab.full-tab[arm]).to_numpy();b=d[ix].mean(1)
                diffrows.append(dict(subset=name,endpoint=metric,control=arm,n_proteins=n,mean_delta=float(d.mean()),median_delta=float(np.median(d)),protein_paired_bootstrap_low=float(np.quantile(b,.025)),protein_paired_bootstrap_high=float(np.quantile(b,.975)),n_full_greater=int((d>0).sum()),n_equal=int((d==0).sum()),n_full_lower=int((d<0).sum()),bootstrap_seed=20261004,n_bootstrap=10000,causal_phi=False))
    pd.DataFrame(summary).to_csv(ROOT/'results'/'WIDE_GROUP_SUMMARY.csv',index=False);pd.DataFrame(diffrows).to_csv(ROOT/'results'/'WIDE_PAIRED_PROTEIN_BOOTSTRAP.csv',index=False)
    print(pd.DataFrame(summary).query("subset=='new25'").to_string(index=False));print(pd.DataFrame(diffrows).query("subset=='new25' and control in ['sequence','frozen_flat','frozen_perm']").to_string(index=False))
    (ROOT/'results'/'WIDE_QA.json').write_text(json.dumps(dict(n_proteins=28,n_new_proteins=25,n_exposed=3,n_tasks=42,conditions=ARMS,n_contexts=16,n_words=1024,n_protein_bootstrap=10000,universal_base_count_predictors=['A','C','G','T'],score_mapping='observed7mer occurrence mean',word_pairs_not_bio_n=True,not_phi_causal_comparison=True,model_and_input_checkpoints_verified=True),indent=2),encoding='utf-8')

if __name__=='__main__':run()
