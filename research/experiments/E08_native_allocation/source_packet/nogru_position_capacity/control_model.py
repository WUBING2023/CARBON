import torch
from MARCP.models.marcp_model import MARCP

class NoGRUPositionCapacity(MARCP):
    def __init__(self):
        super().__init__(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True)
        self.mixer.position_mlp=torch.nn.Sequential(torch.nn.Conv1d(130,216,1),torch.nn.GELU(),torch.nn.Conv1d(216,128,1))
        torch.nn.init.zeros_(self.mixer.position_mlp[2].weight)
        torch.nn.init.zeros_(self.mixer.position_mlp[2].bias)
        self.position_profile_bias=torch.nn.Parameter(torch.zeros(101))
        t=(torch.arange(101,dtype=torch.float32)-50)/50
        self.register_buffer('position_coordinates',torch.stack([t,t.square()]).unsqueeze(0),persistent=False)
    def _backbone(self,x):
        h=self.mixer.residual_proj(self.motif(x.transpose(1,2)))
        coord=self.position_coordinates.to(dtype=h.dtype).expand(h.shape[0],-1,-1)
        h=h+self.mixer.position_mlp(torch.cat([h,coord],dim=1))
        return self.mixer.dropout(self.mixer.norm(h.transpose(1,2))).transpose(1,2)
    def forward(self,x):
        result=super().forward(x)
        result['logits_profile']=result['logits_profile']+self.position_profile_bias.unsqueeze(0)
        return result
