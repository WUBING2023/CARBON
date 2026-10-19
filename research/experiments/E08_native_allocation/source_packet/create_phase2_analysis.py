from pathlib import Path
R=Path(__file__).resolve().parent;O=R/'phase2_ambiguity'
s=(R/'analyze_context.py').read_text(encoding='utf-8').split('df=pd.DataFrame(rows);')[0]
s=s.replace("SCENES=['all_core','identical7','identical7_phi05']","SCENES=['all_core','ambiguous_S','ambiguous_S_phi05','confident_S']")
s=s.replace('default_rng(20261004)','default_rng(20261005)')
s=s.replace("spans=sorted(set((v['start'],v['end']) for v in history))","spans=sorted(set((v['chromosome'],v['start'],v['end']) for v in history))")
s=s.replace("old_names=set(v['gene'] for v in history if v.get('gene'))","old_names=set(json.loads((ROOT/'HISTORICAL_EXPOSURE.json').read_text(encoding='utf-8'))['gene_names'])")
s=s.replace("d=ROOT/'remote/data'/task","d=ROOT/'remote/data'/task/spec['chromosome']")
s=s.replace("exposed_locus=any(v['start']<end and v['end']>start for start,end in spans)","exposed_locus=any(v['chromosome']==ch and v['start']<end and v['end']>start for ch,start,end in spans)")
s=s.replace("'start':v['start']","'chromosome':v['chromosome'],'start':v['start']")
s+='\n'+(O/'analyze_tail.py').read_text(encoding='utf-8')
(O/'analyze_context.py').write_text(s,encoding='utf-8');compile(s,'analyze_context.py','exec');print('Created phase2 analysis')
