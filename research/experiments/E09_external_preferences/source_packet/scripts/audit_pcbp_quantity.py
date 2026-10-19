import h5py,numpy as np,json
from pathlib import Path
root=Path('/data/CARBON_PROJECT/数据/h5data');result=[]
for task in ['PCBP1_HepG2_ENCSR256CHX','PCBP1_K562_ENCSR922WJV']:
    total=0;eligible=0
    for chrom in [8,15]:
        p=root/task/'_5ch'/f'{task}_chr{chrom}_5ch.h5'
        with h5py.File(p) as h:x=h['X'][:]
        good=np.isfinite(x).all((1,2))&np.isin(x[:,:,:4],[0,1]).all((1,2))&(x[:,:,:4].sum(2)==1).all(1);total+=len(x);eligible+=int(good.sum())
    result.append(dict(task=task,n_total=total,n_qualified=eligible,enough_fixed16_contexts=eligible>=16))
out=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/02_external_preferences/PCBP1_SOURCE_QUANTITY.json');out.write_text(json.dumps(result,indent=2));print(json.dumps(result))
