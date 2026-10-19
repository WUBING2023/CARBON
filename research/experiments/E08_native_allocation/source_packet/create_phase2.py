from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'phase2_ambiguity'
s=(ROOT/'remote_context.py').read_text(encoding='utf-8')
s=s.replace('New chr1 inference','New chrX/Y inference')
s=s.replace("OUT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context')","OUT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/phase2_ambiguity')")
s=s.replace("  if row[0]=='chr2' and genes:break\n",'')
s=s.replace("if row[0]=='chr1' and row[2]=='gene':","if row[0] in ['chrX','chrY'] and row[2]=='gene':")
s=s.replace("'start':int(row[3])-1","'chromosome':row[0],'start':int(row[3])-1")
s=s.replace('GENCODE v29 chr1','GENCODE v29 chrX/Y').replace('gene_annotation_chr1_sha256','gene_annotation_chrXY_sha256').replace("'genes_chr1.json'","'genes_chrXY.json'")
start=s.index("for task in M['tasks']:")
end=s.index("source['total_elapsed_seconds']")
block=s[start:end]
lines=block.splitlines()
block=lines[0]+"\n for chrom in ['X','Y']:\n"+'\n'.join(' '+l for l in lines[1:])+'\n'
block=block.replace("task+'_chr1_5ch.h5'","task+'_chr'+chrom+'_5ch.h5'").replace("task+'_chr1.h5'","task+'_chr'+chrom+'.h5'")
block=block.replace("if v['strand']==str(row['strand'])","if v['chromosome']=='chr'+chrom and v['strand']==str(row['strand'])")
block=block.replace("task+'|1|'","task+'|'+chrom+'|'")
block=block.replace('permutations(task,1,len(x))','permutations(task,chrom,len(x))')
block=block.replace("od=OUT/'data'/task;","od=OUT/'data'/task/chrom;")
block=block.replace("{'task':task,'motif':word", "{'task':task,'chromosome':chrom,'motif':word")
block=block.replace("task=task,", "task=task,chromosome=chrom,")
block=block.replace("assert len(df)==len(x);assert x.shape==(len(df),101,5)","assert len(df)==len(x);assert x.shape==(len(df),101,5)\n  assert np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(ctl).all() and (y>=0).all() and (ctl>=0).all()")
block=block.replace("checkpoint_sources=[]", "checkpoint_sources=[]\n  transform_sources={a:hashlib.sha256(np.ascontiguousarray(transformed(z['x'],pm,a)).tobytes()).hexdigest() for a in M['arms']}")
block=block.replace("'source_windows':len(x),'sequence_eligible'", "'transform_sha256':transform_sources,'permutation_sha256':hashlib.sha256(pm.tobytes()).hexdigest(),'source_windows':len(x),'sequence_eligible'")
s=s[:start]+block+s[end:]
(OUT/'remote_context.py').write_text(s,encoding='utf-8')
r=(ROOT/'run_remote.py').read_text(encoding='utf-8').replace("REMOTE='/data/CARBON_PROJECT/supervised_expansion_20261004/01_context'","REMOTE='/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/phase2_ambiguity'")
(OUT/'run_remote.py').write_text(r,encoding='utf-8')
d=(ROOT/'download_remote.py').read_text(encoding='utf-8').replace("/supervised_expansion_20261004/01_context'","/supervised_expansion_20261004/01_context/phase2_ambiguity'")
(OUT/'download_remote.py').write_text(d,encoding='utf-8')
compile(s,str(OUT/'remote_context.py'),'exec')
print('Created phase2 inference source; syntax valid')
