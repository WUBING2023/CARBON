from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'nogru_matched';O=R/'nogru_cross_protein'
src=(S/'remote_train.py').read_text().replace('/01_context/nogru_matched','/01_context/nogru_cross_protein').replace("['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']","['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']").replace('max_batch_seconds\':5400','max_batch_seconds\':2100').replace('>5400','>2100').replace('Registered90min','Registered35min')
(O/'remote_train.py').write_text(src,encoding='utf-8')
(O/'run_remote.py').write_text((S/'run_remote.py').read_text().replace('/01_context/nogru_matched','/01_context/nogru_cross_protein'),encoding='utf-8')
