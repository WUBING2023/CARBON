from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'nogru_matched';O=R/'nogru_position_bias';C=R/'nogru_cross_protein'
cls='''
class NoGRUPosition(MARCP):
    def __init__(self):
        super().__init__(in_channels=5,dropOut=0.2,control=False,ablation_no_gru=True)
        self.position_profile_bias=torch.nn.Parameter(torch.zeros(101))
    def forward(self,x):
        result=super().forward(x)
        result['logits_profile']=result['logits_profile']+self.position_profile_bias.unsqueeze(0)
        return result
'''
src=(S/'remote_train.py').read_text().replace('/01_context/nogru_matched','/01_context/nogru_position_bias').replace("'max_batch_seconds':5400","'max_batch_seconds':2100").replace('>5400','>2100').replace('Registered90min','Registered35min')
src=src.replace('    from MARCP.models.marcp_model import MARCP\n','    from MARCP.models.marcp_model import MARCP\n'+''.join('    '+v+'\n' for v in cls.strip().splitlines()))
src=src.replace('model=MARCP(in_channels=5,dropOut=0.2,control=False,ablation_no_gru=True)','model=NoGRUPosition()').replace('shared={k:full.state_dict()[k].clone() for k in model.state_dict()}','shared={k:(full.state_dict()[k].clone() if k in full.state_dict() else torch.zeros_like(v)) for k,v in model.state_dict().items()}').replace('native noGRU: projection+GRU removed; residual/SE/heads retained','native noGRU+101 zero-init learned absolute-position profile logits bias; count/SE/heads retained')
(O/'remote_train.py').write_text(src,encoding='utf-8')
for filename in ['run_remote.py','download_remote.py','run_infer.py']:(O/filename).write_text((S/filename).read_text().replace('/01_context/nogru_matched','/01_context/nogru_position_bias'),encoding='utf-8')
inf=(S/'remote_infer.py').read_text().replace("O=P/'nogru_matched'","O=P/'nogru_position_bias'").replace("os.environ['PYTHONDONTWRITEBYTECODE']='1'","os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''")
inf=inf.replace('from MARCP.models.marcp_model import MARCP','from MARCP.models.marcp_model import MARCP\n'+cls)
inf=inf.replace('model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True)','model=NoGRUPosition()').replace('{k:full.state_dict()[k].clone() for k in model.state_dict()}','{k:(full.state_dict()[k].clone() if k in full.state_dict() else torch.zeros_like(v)) for k,v in model.state_dict().items()}').replace("z['nogru_seed'","z['nogru_pos_seed'")
(O/'remote_infer.py').write_text(inf,encoding='utf-8')
for filename in ['remote_infer.py','run_infer.py','download_remote.py']:
 cs=(S/filename).read_text().replace('/01_context/nogru_matched','/01_context/nogru_cross_protein').replace("O=P/'nogru_matched'","O=P/'nogru_cross_protein'").replace("['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']","['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']")
 if filename=='remote_infer.py':cs=cs.replace("os.environ['PYTHONDONTWRITEBYTECODE']='1'","os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''")
 (C/filename).write_text(cs,encoding='utf-8')
print('Created registered zero-position-bias training and strict inference; HNR inference is separate.')
