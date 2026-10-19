from pathlib import Path
R=Path(__file__).resolve().parents[1]
source=(R/'scripts'/'nogru_external.py').read_text()
source=source.replace("'nogru_matched'/'outputs'","'nogru_cross_protein'/'outputs'").replace("['PTBP1_HepG2_ENCSR384KAN','PTBP1_K562_ENCSR981WKN']","['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU']").replace("'nogru_external_outputs'","'hnrnpk_nogru_external_outputs'")
(R/'scripts'/'hnrnpk_nogru_external.py').write_text(source,encoding='utf-8')
print('Prepared unchanged HNRNPK controlled-query method')
