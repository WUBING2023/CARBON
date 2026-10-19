from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'nogru_matched';O=R/'nogru_cross_protein'
s=(S/'analyze_context.py').read_text().replace('PTBP1','HNRNPK')
(O/'analyze_context.py').write_text(s,encoding='utf-8')
