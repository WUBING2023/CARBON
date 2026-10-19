from pathlib import Path
R=Path(__file__).resolve().parents[1]
src=(R/'scripts'/'cross_audit_context.py').read_text().replace("NG=A/'nogru_matched'","NG=A/'nogru_cross_protein'").replace("['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']","['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']").replace("group=='PTBP1'","group=='HNRNPK'").replace('all_core_adjusted_PTBP1.npz','all_core_adjusted_HNRNPK.npz').replace('CROSS_AUDIT_01_','CROSS_AUDIT_01_HNR_')
first=src.index("bs=np.load(NG/'bootstrap'")
last=src.index('pd.DataFrame(checks)',first)
src=src[:first]+'''bs=np.load(NG/'bootstrap'/'all_core_adjusted_HNRNPK.npz');means=bs['means'];methods=list(bs['methods']);delta=means[:,methods.index('aligned_ensemble')]-means[:,methods.index('nogru_ensemble')]
contrast=pd.read_csv(NG/'CONTRASTS.csv');row=contrast.query("scene=='all_core' and cohort=='adjusted' and group=='HNRNPK' and seed=='ensemble' and control=='nogru'").iloc[0]
assert abs(float(row['gene_equal_delta'])-(float(p['aligned_ensemble_loss'])-float(p['nogru_ensemble_loss'])))<1e-11
assert np.allclose(np.quantile(delta,[.025,.975]),row[['bonf1_low','bonf1_high']].to_numpy(float),atol=1e-11)
savedDelta=delta.copy()
''' +src[last:]
src=src.replace('bootstrap_10000_max_abs_error=float(abs(delta-savedDelta).max())','CI_from_saved_arm_draws_independently_recomputed=True,bootstrap_draw_generation_repeated=False')
(R/'scripts'/'cross_audit_hnr_context.py').write_text(src,encoding='utf-8')
src=(R/'scripts'/'audit_nogru_initialization.py').read_text().replace('/01_context/nogru_matched','/01_context/nogru_cross_protein').replace("['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']","['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']").replace('CROSS_AUDIT_01_INITIALIZATION.json','CROSS_AUDIT_01_HNR_INITIALIZATION.json')
(R/'scripts'/'audit_hnr_initialization.py').write_text(src,encoding='utf-8')
print('Prepared independent HNR primary/seed/cell and initialization audits')
