from pathlib import Path
R=Path(__file__).resolve().parent;O=R/'nogru_position_capacity';S=R/'nogru_cross_protein'
s=(S/'remote_pilot.py').read_text().replace('/01_context/nogru_cross_protein','/01_context/nogru_position_capacity').replace("['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID','RBFOX2_K562_ENCSR756CKJ','RBFOX2_HepG2_ENCSR987FTF']","['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']").replace('from MARCP.models.marcp_model import MARCP','from MARCP.models.marcp_model import MARCP\nsys.path.insert(0,str(O));from control_model import NoGRUPositionCapacity').replace('model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True)','model=NoGRUPositionCapacity()').replace('shared={k:full.state_dict()[k].clone() for k in model.state_dict()}','shared={k:(full.state_dict()[k].clone() if k in full.state_dict() else v.clone()) for k,v in model.state_dict().items()}').replace("'PILOT_PROTOCOL.md'","'PROTOCOL.md'")
s=s.replace('params=sum(p.numel() for p in model.parameters());h=state_hash(model.state_dict());model.to(\'cuda:0\')',"""params=sum(p.numel() for p in model.parameters());h=state_hash(model.state_dict());assert params==100871
 base=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);base.load_state_dict({k:model.state_dict()[k].clone() for k in base.state_dict()},strict=True);base.eval();model.eval()
 qa=torch.linspace(-1,1,16*101*5).reshape(16,101,5)
 with torch.no_grad():
  aa=base(qa);bb=model(qa);initial_qa={k:float((aa[k]-bb[k]).abs().max()) for k in aa}
 assert all(v==0 for v in initial_qa.values());del base;model.to('cuda:0')""")
s=s.replace("'parameters':params,","'parameters':params,'initial_output_qa':initial_qa,")
s=s.replace('upper12runs100epochs_seconds','upper6runs100epochs_seconds').replace("'feasible_under60min'","'feasible_under45min'").replace('<=3600','<=2700')
(O/'remote_pilot.py').write_text(s,encoding='utf-8')
tr=(O/'run_remote.py').read_text().replace("'remote_train.py'","'remote_pilot.py'").replace('r/"remote_train.py"','r/"remote_pilot.py"').replace('REMOTE_LOG.jsonl','PILOT_REMOTE_LOG.jsonl').replace('TRANSPORT_RUN.json','TRANSPORT_PILOT.json')
(O/'run_pilot.py').write_text(tr,encoding='utf-8')
print('Capacity pilot freezes final architecture, identical initial outputs and all100budget costs.')
