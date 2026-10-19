from pathlib import Path
R=Path(__file__).resolve().parent
s=(R/'verify_weighted_refits.py').read_text('utf8').split('# Gene x RBP')[0]
s=s.replace("R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv'","R/'STRICT396_WITH_ANNOTATIONS.csv'")
s=s.replace("[['task','strand']]","[['task','strand','NCBI_annotation_region']]")
s=s.replace('20261014','20261016').replace('WEIGHTED_REFIT_SENSITIVITY_200','WEIGHTED_REFIT_REGION_SENSITIVITY_200').replace('WEIGHTED_REFIT_SAMPLES_200','WEIGHTED_REFIT_REGION_SAMPLES_200')
(R/'verify_region_sensitivity.py').write_text(s,'utf8')
(R/'REGION_SENSITIVITY_PLAN.md').write_text('Secondary exploratory sensitivity added after primary descriptive results because exact NCBI GRCh38 consequences became available. Add NCBI region category to the frozen nuisance controls; retain all 396 events including unavailable category. Coding/splice/UTR/noncoding exon/intronic/other/unavailable precedence is declared in finalize_annotations.py. Refit the same centerPhi slopes for the same eight outcomes, using 200 positive exponential two-way gene x RBP weighted refits, seed20261016. Do not alter any primary endpoint, plot or cohort based on this extension.','utf8')
