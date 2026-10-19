from pathlib import Path
R=Path(__file__).resolve().parents[1];s=(R/'scripts'/'cross_audit_context.py').read_text()
s=s.replace("NG=A/'nogru_matched'","NG=A/'nogru_position_capacity'").replace("'permuted','nogru']","'permuted','nogru','nogru_position_capacity']")
s=s.replace("profiles={}\n        for arm", "oldng=np.load(A/'nogru_matched'/'remote'/'data'/task/chrom/'data.npz')\n        assert np.array_equal(oldng['index'],z['index'])\n        profiles={}\n        for arm")
s=s.replace("store=zz if arm=='nogru' else z","store=zz if arm=='nogru_position_capacity' else (oldng if arm=='nogru' else z)")
s=s.replace("['aligned_ensemble_loss','nogru_ensemble_loss']","['aligned_ensemble_loss','nogru_position_capacity_ensemble_loss']").replace("methods.index('nogru_ensemble')","methods.index('nogru_position_capacity_ensemble')")
s=s.replace('CROSS_AUDIT_01_','CROSS_AUDIT_01_PC_')
(R/'scripts'/'cross_audit_capacity_context.py').write_text(s,encoding='utf-8')
print('Prepared independent updated-control native readback')
