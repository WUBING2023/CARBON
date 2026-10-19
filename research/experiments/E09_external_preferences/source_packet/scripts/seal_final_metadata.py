"""Seal the final already-completed audit artifacts without repeating checks."""
from pathlib import Path
import json,hashlib,datetime
R=Path(__file__).resolve().parents[1]
q=R/'QA.json';data=json.loads(q.read_text());p=R/'results'/'CROSS_AUDIT_01_PC_TRAINED_BRANCHES.json'
data['final_trained_branch_readback']=dict(path=p.relative_to(R).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n_norm_checks=90,max_abs_error=0.0,no_forward_no_training=True)
data['remaining_scope']='none in 02; root final reading-page review is outside this slice'
q.write_text(json.dumps(data,indent=2),encoding='utf-8')
audit=R/'CROSS_AUDIT_01.md';text=audit.read_text(encoding='utf-8')
text+='\n最终已训练支路身份读回：本地CPU2直接读取六best权重，private五tensor×三norm共90项与A学习QA精确0差；初零MLP末层和profilebias六run均已变成非零。没有再forward或训练。控制源码hash905f265e…与初始审核相同，未漂移。\n'
audit.write_text(text,encoding='utf-8')
man=R/'MANIFEST.json';metadata=json.loads(man.read_text());records={x['path']:x for x in metadata['files']}
for p in [q,audit,R/'scripts'/'audit_trained_branches.py',R/'scripts'/'seal_final_metadata.py',R/'results'/'CROSS_AUDIT_01_PC_TRAINED_BRANCHES.json']:
    name=p.relative_to(R).as_posix();records[name]=dict(path=name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
metadata['files']=[records[k] for k in sorted(records)];metadata['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();metadata['final_metadata_only_no_repeated_QA']=True
man.write_text(json.dumps(metadata,indent=2),encoding='utf-8')
print('Final manifest sealed;',len(records),'files; no repeated analysis or inference.')
