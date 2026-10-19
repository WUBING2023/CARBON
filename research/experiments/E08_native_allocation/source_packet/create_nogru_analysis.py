from pathlib import Path
R=Path(__file__).resolve().parent;O=R/'nogru_matched'
s=(R/'phase2_ambiguity/analyze_context.py').read_text(encoding='utf-8').split('df=pd.DataFrame(rows)')[0]
s=s.replace("ARMS=['aligned','sequence_only','mean','permuted']","ARMS=['aligned','sequence_only','mean','permuted','nogru']")
s=s.replace("SCENES=['all_core','ambiguous_S','ambiguous_S_phi05','confident_S']","SCENES=['all_core','long_gap','identical_local_input']")
s=s.replace('default_rng(20261005)','default_rng(20261006)')
s=s.replace("ROOT/'HISTORICAL_EXPOSURE.json'","ROOT.parent/'HISTORY_REPLAY_ALL_CHROMOSOMES.json'")
s=s.replace("source=json.loads((ROOT/'remote/SOURCE_MANIFEST.json').read_text(encoding='utf-8'))","source=json.loads((ROOT/'remote/INFERENCE_MANIFEST.json').read_text(encoding='utf-8'))")
s=s.replace("d=ROOT/'remote/data'/task/spec['chromosome']\n with np.load(d/'data.npz') as f:z={k:f[k] for k in f.files}\n meta=json.loads((d/'metadata.json').read_text(encoding='utf-8'))", "d=ROOT.parent/'remote/data'/task if spec['chromosome']=='1' else ROOT.parent/'phase2_ambiguity/remote/data'/task/spec['chromosome']\n with np.load(d/'data.npz') as f:z={k:f[k] for k in f.files}\n meta=json.loads((d/'metadata.json').read_text(encoding='utf-8'))\n with np.load(ROOT/'remote/data'/task/spec['chromosome']/'data.npz') as f:nz={k:f[k] for k in f.files}\n assert np.array_equal(nz['index'],z['index']) and np.array_equal(nz['centers'],z['centers'])\n z.update({k:v for k,v in nz.items() if k.startswith('nogru_')})")
s=s.replace("row={'task':task", "identical_local=8<=c1<=92 and 8<=c2<=92 and np.array_equal(z['x'][i,c1-8:c1+9,:4],z['x'][i,c2-8:c2+9,:4]) and np.max(np.abs(z['x'][i,c1-8:c1+9,4]-z['x'][i,c2-8:c2+9,4]))<=1e-6\n  row={'range':('chr1' if spec['chromosome']=='1' else 'XY'),'long_gap':c2-c1>=31,'identical_local_input':identical_local,'task':task")
tail=(R/'phase2_ambiguity/analyze_tail.py').read_text(encoding='utf-8')
tail=tail.replace("for co in ['raw','adjusted']:","for co in ['raw','adjusted']:")
tail=tail.replace("for scene in SCENES:\n for cohort", "for scene in SCENES:\n for cohort")
tail=tail.replace("for group,tt in [(t,[t]) for t in tasks]+[(p,[t for t in tasks if t.startswith(p+'_')]) for p in ['PTBP1','HNRNPK','RBFOX2']]:", "for group,tt,ranges in [(t,[t],['chr1','XY']) for t in tasks]+[(g,tasks,rs) for g,rs in [('PTBP1',['chr1','XY']),('PTBP1_chr1',['chr1']),('PTBP1_XY',['XY'])]]:")
tail=tail.replace("sub=sd[sd['task'].isin(tt)];gg=g[g.index.get_level_values('task').isin(tt)];names=", "sub=sd[sd['task'].isin(tt)&sd['range'].isin(ranges)];gg=sub.groupby(['task','gene_id'],sort=True)[cols].mean();names=")
tail=tail.replace("['sequence_only','mean','permuted']", "['nogru','sequence_only','mean','permuted']")
tail=tail.replace('for fam in [1,3,27,54]:','for fam in [1,12,36]:')
tail=tail.replace("'rng_seed':20261005", "'rng_seed':20261006")
tail=tail.replace("'primary':'PTBP1 / ambiguous_S / adjusted; protein-level same geneID clusters jointly across cells;3controls'", "'primary':'PTBP1 combinedchr1XY / allcore / adjusted, Full−NoGRU;95%CI, samegeneIDclusters joint acrosscells'")
tail=tail.replace("'supporting':'protein×scene×control27; task×scene×control54. confident_S complement descriptive only, not confirmatory family.'", "'supporting':'3scenes×4controls Bonf12 for combinedrange; Bonf36 conservativeall3ranges. Bothranges also reportedseparately. ModelsseedsnotRNAreplicates.'")
tail=tail.replace("(c['group']=='PTBP1')", "(c['group'].str.startswith('PTBP1'))")
s+='\n'+tail
(O/'analyze_context.py').write_text(s,encoding='utf-8');compile(s,'analyze_context.py','exec');print('Created NoGRU analysis source')
