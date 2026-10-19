from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'nogru_matched';O=R/'nogru_position_capacity'
s=(S/'analyze_context.py').read_text()
s=s.replace("'permuted','nogru'","'permuted','nogru','nogru_position_capacity'").replace('RNG=np.random.default_rng(20261006)','RNG=np.random.default_rng(20261008)')
s=s.replace("z.update({k:v for k,v in nz.items() if k.startswith('nogru_')})","""z.update({k:v for k,v in nz.items() if k.startswith('nogru_position_capacity_')})
 with np.load(ROOT.parent/'nogru_matched/remote/data'/task/spec['chromosome']/'data.npz') as f:oz={k:f[k] for k in f.files}
 assert np.array_equal(oz['index'],z['index']) and np.array_equal(oz['centers'],z['centers'])
 z.update({k:v for k,v in oz.items() if k.startswith('nogru_seed')})""")
s=s.replace("for control in ['nogru','sequence_only','mean','permuted']:","for control in ['nogru_position_capacity','nogru','sequence_only','mean','permuted']:")
s=s.replace('for fam in [1,12,36]:','for fam in [1,2,15,45]:')
s=s.replace("'rng_seed':20261006","'rng_seed':20261008").replace('Full−NoGRU;95%CI','Full−NoGRU_position_capacity;95%CI').replace('3scenes×4controls Bonf12','3scenes×5controls Bonf15').replace('Bonf36','Bonf45')
s=s.replace("SCENES=['all_core','long_gap','identical_local_input']","SCENES=['all_core','long_gap','identical_local_input','distance_equal','distance_le2','distance_le5']")
s=s.replace("df['ambiguous_S']=", "df['distance_difference']=np.abs(np.abs(df.c1-50)-np.abs(df.c2-50))\ndf['distance_equal']=df.distance_difference==0\ndf['distance_le2']=df.distance_difference<=2\ndf['distance_le5']=df.distance_difference<=5\ndf['ambiguous_S']=")
# Save the prespecified second diagnostic contrast without altering any model losses.
s=s.replace('contrasts.append(rec)',"""contrasts.append(rec)
   for seed in ['ensemble','seed11','seed23','seed47']:
    ai=METHODS.index('nogru_position_capacity_'+seed);ci=METHODS.index('nogru_'+seed);vv=bs[:,ai]-bs[:,ci];vv=vv[np.isfinite(vv)]
    rec={'scene':scene,'cohort':cohort,'group':group,'seed':seed,'control':'nogru','comparison':'position_capacity_minus_nogru','gene_equal_delta':point[ai]-point[ci],'pair_pooled_delta':(sub['nogru_position_capacity_'+seed+'_'+co+'_loss']-sub['nogru_'+seed+'_'+co+'_loss']).mean(),'genes':len(names),'loci':len(sub),'valid_bootstraps':len(vv)}
    for fam in [1,2,15,45]:
     rec['bonf'+str(fam)+'_low']=np.quantile(vv,.05/(2*fam)) if len(vv) else None;rec['bonf'+str(fam)+'_high']=np.quantile(vv,1-.05/(2*fam)) if len(vv) else None
    contrasts.append(rec)""")
(O/'analyze_context.py').write_text(s,encoding='utf-8')
print('Registered capacity-control losses; seed20261008 and full paired contrasts preserved.')
