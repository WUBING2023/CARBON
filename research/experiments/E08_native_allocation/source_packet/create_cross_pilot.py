from pathlib import Path
R=Path(__file__).resolve().parent;O=R/'nogru_cross_protein';O.mkdir(exist_ok=True)
s=(R/'nogru_pilot/remote_pilot.py').read_text(encoding='utf-8')
s=s.replace("O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_pilot')","O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_cross_protein')")
s=s.replace("['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']","['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID','RBFOX2_K562_ENCSR756CKJ','RBFOX2_HepG2_ENCSR987FTF']")
s=s.replace("sha(O/'PROTOCOL.md')","sha(O/'PILOT_PROTOCOL.md')").replace('upper6runs100epochs_seconds','upper12runs100epochs_seconds').replace('feasible_under90min','feasible_under60min').replace('<=5400','<=3600')
(O/'remote_pilot.py').write_text(s,encoding='utf-8')
s=(R/'run_remote.py').read_text(encoding='utf-8').replace('/01_context\'','/01_context/nogru_cross_protein\'').replace('remote_context.py','remote_pilot.py').replace('PROTOCOL.md','PILOT_PROTOCOL.md').replace('REMOTE_LOG.jsonl','PILOT_LOG.jsonl').replace('TRANSPORT_RUN.json','TRANSPORT_PILOT.json')
(O/'run_pilot.py').write_text(s,encoding='utf-8');print('Cross-protein pilot sources created')
