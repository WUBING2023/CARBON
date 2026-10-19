"""只读核对已保存证据，写出本次检查结果；不训练/推理，不覆盖实验表。

检查内容：原指标→消融配对表；逐pair重新计算target/loss及gene/cell聚合；
RNAcompete全部5mer的Spearman/同组成排序；疾病完整队列及权重位置。
这证明已检查的来源链/算术一致，不是排除一切科研不端的法证证明。
"""
from pathlib import Path
import csv,json,math,itertools,collections,hashlib,sys
import numpy as np
from scipy.stats import spearmanr
OUT=Path(__file__).resolve().parents[1]
EXP=OUT/'research/experiments';SERVER=EXP/'E03_ablation/original_records/全量消融运行'
import argparse
cli=argparse.ArgumentParser(description=__doc__);cli.add_argument('--output',type=Path,default=OUT/'outputs/frozen_results');args=cli.parse_args();AUDIT=args.output;AUDIT.mkdir(parents=True,exist_ok=True)
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if str(p).endswith('.tsv') else ','))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def avg(v):return math.fsum(v)/len(v)
def metrics(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return {x['Metric']:float(x['Value']) for x in csv.DictReader(f,delimiter='\t')}
def write_csv(p,data,fields=None):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(data[0]));w.writeheader();w.writerows(data)

def ablations():
    data=rows(EXP/'E03_ablation/paper_data/Supplementary_Data_S2_ablation_paired_deltas.csv')
    assert len(data)==1110 and len({r['task'] for r in data})==222
    keys=['Mean_PCC','Spearman','Signal_AUC','Signal_AP','Peak_AUC','Peak_AP'];maxerror=0.;checks=0;runindex={}
    index=rows(OUT/'research/catalog/ABLATION_RUN_INDEX.csv');assert len(index)==1332
    runindex={(r['task'],r['arm']):r for r in index}
    for r in data:
        for arm,prefix in [(r['ablation'],''),('baseline','baseline_')]:
            run=runindex[(r['task'],arm)];original=metrics(OUT/run['packaged_metric'])
            for k in keys:
                error=abs(original[k]-float(r[prefix+k]));maxerror=max(maxerror,error);assert error<2e-12,(r['task'],arm,k,error);checks+=1
    summaries=[]
    for arm in sorted({r['ablation'] for r in data}):
        a=[r for r in data if r['ablation']==arm]
        summaries.append(dict(arm=arm,n_tasks=len(a),mean_delta_spearman=avg([float(r['Spearman'])-float(r['baseline_Spearman']) for r in a]),spearman_decreased=sum(float(r['Spearman'])<float(r['baseline_Spearman']) for r in a),mean_delta_pcc=avg([float(r['Mean_PCC'])-float(r['baseline_Mean_PCC']) for r in a]),pcc_decreased=sum(float(r['Mean_PCC'])<float(r['baseline_Mean_PCC']) for r in a)))
    write_csv(AUDIT/'ABLATION_RECOMPUTED_SUMMARY.csv',summaries)
    return dict(rows=len(data),tasks=222,paired_metric_comparisons=checks,maximum_absolute_error=maxerror,runs=len(runindex),summaries=summaries,scope='Saved original per-run metrics versus paper paired table; not rerunning training or metric extraction from raw signal.')

def benchmark():
    data=rows(OUT/'research/catalog/BENCHMARK_WEIGHT_AVAILABILITY.csv')
    counts={m:sum(r[m+'_exists']=='True' for r in data) for m in ['CARBON','iDeepB','RBPNet']}
    assert counts=={'CARBON':221,'iDeepB':219,'RBPNet':222},counts
    return dict(tasks=len(data),weight_counts=counts,scope='Canonical archived weight identities, not bundled files or a fair benchmark rerun.')

def allocation():
    p=EXP/'E08_native_allocation/source_packet/nogru_position_capacity'
    data=rows(p/'PAIR_RESULTS_ALL.csv');primary=[r for r in data if r['rbp']=='PTBP1' and r['all_core']=='True' and r['unique_gene']=='True' and r['adj_good']=='True']
    assert len(primary)==1851
    arms=['aligned_ensemble','nogru_ensemble','nogru_position_capacity_ensemble'];group=collections.defaultdict(list);max_q=0.;max_loss=0.
    for r in primary:
        scale=float(r['signal_total'])/float(r['input_total']);a=float(r['y1'])-scale*float(r['input1']);b=float(r['y2'])-scale*float(r['input2']);assert a>=-1e-10 and b>=-1e-10 and a+b>0
        q=a/(a+b);max_q=max(max_q,abs(q-float(r['q_adjusted'])))
        for arm in arms:
            share=float(r[arm+'_mass1'])/(float(r[arm+'_mass1'])+float(r[arm+'_mass2']));loss=(share-q)**2
            max_loss=max(max_loss,abs(loss-float(r[arm+'_adjusted_loss'])));group[(arm,r['cell'],r['gene_id'])].append(loss)
    reference={r['method']:float(r['gene_equal_mse']) for r in rows(p/'METHODS.csv') if r['scene']=='all_core' and r['cohort']=='adjusted' and r['group']=='PTBP1'}
    point={}
    for arm in arms:
        cells=collections.defaultdict(list)
        for (method,cell,gene),values in group.items():
            if method==arm:cells[cell].append(avg(values))
        point[arm]=avg([avg(v) for v in cells.values()]);assert abs(point[arm]-reference[arm])<3e-10,(arm,point[arm],reference[arm])
    assert max_q<3e-10 and max_loss<3e-10
    return dict(pairs=len(primary),genes=len({r['gene_id'] for r in primary}),cell_gene_units=len({(r['cell'],r['gene_id']) for r in primary}),maximum_q_error=max_q,maximum_loss_error=max_loss,cell_equal_gene_mse=point,full_vs_local_pos_relative_reduction_percent=100*(1-point[arms[0]]/point[arms[2]]),scope='Reconstructed input-adjusted targets, prediction shares, pair losses and gene/cell aggregation from complete saved pair table; not fresh neural inference or bootstrap CI rerun.')

def external():
    p=EXP/'E09_external_preferences/source_packet/results'
    pred=collections.defaultdict(dict);truth=collections.defaultdict(dict);meta={};kmers=[''.join(x) for x in itertools.product('ACGT',repeat=5)]
    for r in rows(p/'WIDE_ALL_5MER_PREDICTIONS.csv'):
        if r['arm'] not in ['full','sequence','frequency']:continue
        key=(r['task'],r['arm']);pred[key][r['kmer']]=float(r['score']);meta[key]=(r['protein'],r['cell'],r['exposed_protein']=='True')
    for r in rows(p/'WIDE_EXTERNAL_5MER_SCORES.csv'):truth[(r['protein'],r['experiment'])][r['kmer']]=float(r['z5'])
    composition=collections.defaultdict(list)
    for i,k in enumerate(kmers):composition[tuple(k.count(c) for c in 'ACGT')].append(i)
    pairs=np.array([ab for g in composition.values() for ab in itertools.combinations(g,2)],dtype=int);assert len(pairs)==15240
    original={(r['task'],r['arm'],r['experiment']):r for r in rows(p/'WIDE_ALL_METRICS.csv')}
    errors=[];bycell=collections.defaultdict(lambda:[[],[]]);n=0
    for key,values in pred.items():
        protein,cell,exposed=meta[key];pr=np.array([values[k] for k in kmers]);dp=pr[pairs[:,0]]-pr[pairs[:,1]]
        for (name,assay),target in truth.items():
            if name!=protein:continue
            ac=np.array([target[k] for k in kmers]);dt=ac[pairs[:,0]]-ac[pairs[:,1]];valid=dt!=0
            rho=float(spearmanr(pr,ac).statistic);cc=float(np.mean((dp[valid]*dt[valid]>0)+.5*(dp[valid]==0)))
            old=original[(key[0],key[1],assay)];errors.extend([abs(rho-float(old['spearman_all5mer'])),abs(cc-float(old['within_composition_concordance']))]);n+=1
            if not exposed:
                bycell[(protein,cell,key[1])][0].append(rho);bycell[(protein,cell,key[1])][1].append(cc)
    assert max(errors)<1e-12,max(errors)
    byprotein=collections.defaultdict(lambda:[[],[]])
    for (protein,cell,arm),v in bycell.items():
        byprotein[(protein,arm)][0].append(avg(v[0]));byprotein[(protein,arm)][1].append(avg(v[1]))
    result={arm:dict(spearman_all5mer=avg([avg(v[0]) for (protein,a),v in byprotein.items() if a==arm]),within_composition_concordance=avg([avg(v[1]) for (protein,a),v in byprotein.items() if a==arm])) for arm in ['full','sequence','frequency']}
    assert len({protein for protein,arm in byprotein})==25
    return dict(proteins=25,words=1024,same_composition_pairs=15240,task_arm_assay_checks=n,maximum_metric_error=max(errors),recomputed_protein_equal_means=result,full_minus_frequency_spearman=result['full']['spearman_all5mer']-result['frequency']['spearman_all5mer'],full_minus_sequence_concordance_pp=100*(result['full']['within_composition_concordance']-result['sequence']['within_composition_concordance']),scope='All full/sequence/frequency 5mer tables versus external assay values; ties exact, assay then cell then protein aggregation. Does not regenerate RNAcompete7→5 mapping or run neural predictions.')

def disease():
    p=EXP/'E10_natural_variants/disease_annotations';a=rows(p/'DISEASE_ALL13_EVENT_RESULTS.csv');s=json.loads((p/'DISEASE_SUMMARY.json').read_text())
    result=dict(events=len(a),variants=len({r['variant_key'] for r in a}),full_matches=sum(float(r['A_native_credit']) for r in a),sequence_matches=sum(float(r['S_credit']) for r in a),split_counts=dict(collections.Counter(r['split'] for r in a)))
    assert result['events']==s['union_events']==13 and result['variants']==s['union_variants']==8
    assert result['full_matches']==s['union_full_direction_matches']==7 and result['sequence_matches']==s['union_sequence_direction_matches']==10
    result['scope']='Complete saved clinical/tumor subset; preserves unfavorable Full7/13 vs Sequence10/13. Annotation is not independent functional validation.'
    return result

def main():
    result={};failures=[]
    for name,fn in [('ablations',ablations),('benchmark_availability',benchmark),('allocation',allocation),('external_preferences',external),('disease_subset',disease)]:
        try:result[name]=fn();print(name,'PASS',flush=True)
        except Exception as e:result[name]={'status':'FAIL','error':str(e)};failures.append(name);print(name,'FAIL',str(e),flush=True)
    result['check_status']='PASS' if not failures else 'FAIL';result['failed_checks']=failures
    result['not_established']=['No forensic certification that every original measurement is authentic.','No fresh training/inference, no unified three-model benchmark, no new biological validation.','Raw sequencing/reference data not fully content-hashed; server metadata is indexed.']
    dump(AUDIT/'NUMERICAL_RECHECK.json',result)
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
