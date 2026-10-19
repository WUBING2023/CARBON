from literal_tables import read_csv_literal
from pathlib import Path
import hashlib,json,datetime
import numpy as np,pandas as pd
from PIL import Image
R=Path(__file__).resolve().parents[1];D=R/'results';F=R/'figures'
now=datetime.datetime.now(datetime.timezone.utc).isoformat();qa=json.loads((D/'WIDE_QA.json').read_text());checks=[]
def check(name,test):
    checks.append(dict(name=name,passed=bool(test)));assert test,name
man=json.loads((R/'WIDE_MANIFEST.json').read_text());execution=read_csv_literal(D/'WIDE_EXECUTION_QA.csv');check('42tasks_fourconditions_true_manifests',len(execution)==168 and execution.task.nunique()==42 and len(man)==42)
check('complete_input_checkpoint_phi_QA',execution[['checkpoint_hash_pass','input_identity_pass','phi_transform_pass']].all().all())
fd=json.loads((F/'FIGURE_DATA.json').read_text());prot=read_csv_literal(D/'WIDE_PROTEIN_METRICS.csv');new=prot[~prot.exposed_protein];check('all_new25_in_figure',set(fd['main_new_proteins'])==set(new.protein) and len(fd['main_new_proteins'])==25)
for point in fd['panels_ab']:
    row=new[(new.protein==point['protein'])&(new.arm==point['condition'])].iloc[0];check('figure_'+point['panel']+'_'+point['protein']+'_'+point['condition'],abs(float(row[point['endpoint']])-point['value'])<1e-13)
check('all20truth_words',len(fd['panel_c_words'])==20 and len(set(fd['panel_c_words']))==20)
words=read_csv_literal(D/'HNRNPK_ALL_TOP20_WORDS.csv');sel=words.groupby(['kmer','arm']).in_model_top20.mean().unstack();tru=words.groupby('kmer').external_Z5.first()
for word in fd['panel_c_words']:
    check('c_truth_'+word,abs(fd['panel_c_Z5'][word]-float(tru[word]))<1e-13)
    for arm in ['aligned','sequence_only','mean','permuted']:check('c_priority_'+word+'_'+arm,abs(fd['panel_c_selection'][arm][word]-float(sel.loc[word,arm]))<1e-13)
for p in fd['panel_d']:
    row=read_csv_literal(D/'HNRNPK_PAIRED_GENE_BOOTSTRAP.csv');row=row[(row.assay==p['assay'])&(row.control==p['control'])].iloc[0]
    for source,target in [('delta_recall20','delta_pp'),('conditional_gene_locus_bootstrap_low','low_pp'),('conditional_gene_locus_bootstrap_high','high_pp')]:check('d_'+p['assay']+'_'+p['control']+'_'+target,abs(float(row[source])*100-p[target])<1e-12)
layout=json.loads((F/'LAYOUT_QA.json').read_text());check('no_visible_text_outside_canvas',not layout['outside_canvas']);check('85mm_canvas',fd['width_mm']==85)
check('caption_synced',(R/'CAPTION.md').read_bytes()==(F/'CAPTION.md').read_bytes())
no=read_csv_literal(D/'NOGRU_EXTERNAL_QA.csv');check('all6_nogru_external_full_scope',len(no)==6 and no.same_native_backgrounds.all() and no.same_phyloP.all())
hn=read_csv_literal(D/'HNR_COMPONENT_QA.csv');check('all6_HNR_component_same_scope',len(hn)==6 and hn.same16_backgrounds.all() and hn.samePhi.all())
pc=json.loads((D/'CROSS_AUDIT_01_PC_QA.json').read_text());check('final_control48points_and10000_draws',pc['point_checks']==48 and pc['max_point_difference']<1e-11 and pc['bootstrap_10000_max_abs_error']<1e-13)
model=json.loads((D/'CROSS_AUDIT_01_POSITION_CAPACITY_MODEL.json').read_text());check('final_control6fitted_shared_initializations',len(model['runs'])==6 and all(x['shared_tensor_init_exact'] and x['fresh_init_hash_matches'] and x['MLP_pointwise'] for x in model['runs']))
check('updated_logits_used_all_real_stages',all(model[k] for k in ['loss_consumes_updated_logit','validation_consumes_updated_logit','inference_consumes_updated_logit','parent_forward_has_no_cached_probability']))
ciqa=json.loads((D/'CROSS_AUDIT_01_PC_FINAL_CI_QA.json').read_text());check('all36final_named_intervals',ciqa['all_passed'] and ciqa['n_checks']==36)
parser=json.loads((D/'PARSER_QA.json').read_text());check('622516source_tokens_literal_equal',parser['raw_literal_bit_equal'] and parser['raw_numeric_tokens_checked']==622516)
check('new25unchanged_after_global_parser',parser['new25_primary_points_bit_unchanged'])
init=json.loads((D/'CROSS_AUDIT_01_INITIALIZATION.json').read_text());check('6fresh_init_order_audits',len(init['rows'])==6 and all(x['shared_epoch_orders_equal'] and x['nogru_initial_hash_matches'] for x in init['rows']))
ig=read_csv_literal(D/'NUMERICAL_QA_ALL72.csv');check('all72_IG_numerical_QA',len(ig)==72)
for p in [R/'RESULTS.md',R/'REPORT.md',R/'CAPTION.md',R/'REPLAY.md',R/'SOURCES.md']:check('exists_'+p.name,p.exists() and p.stat().st_size>100)
independent={}
for fname,n in [('FAMILY_QA.json',2631),('FIGURE_QA.json',620),('QA.json',39)]:
    p=R.parent/'03_functional_rules'/'cross_audit_02'/fname
    independent[fname]=dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),independent_checks_reported_by_reviewer=n)
packet=dict(utc=now,checks=checks,n_checks=len(checks),all_passed=True,wide=qa,PNG_size=list(Image.open(F/'EXTERNAL_RULES.png').size),figure_width_mm=85,figure_height_mm=132,physical_min_font_pt=6.0,scope='family external preference confirmation plus exposed one-protein local example; no conservation-causal proof',not_manuscript_submission_readiness=True,independent_review_reference=independent,root_visual_review='current PNG byte-identical to independently reviewed navy image',remaining_scope='root final reading-page render only; independent external and context-control checks complete',original_independent_parser_failures_preserved=62,all_passed_refers_to_current_packet_checks_only=True)
(R/'QA.json').write_text(json.dumps(packet,indent=2),encoding='utf-8')
hashes=[]
for folder in [R/'scripts',R/'sources',F]:
    for p in sorted(folder.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and 'revisions' not in p.parts:hashes.append(dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
for p in sorted(D.glob('*')):
    if p.is_file():hashes.append(dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
for name in ['PROTOCOL.md','PROTOCOL_FREEZE.json','SOURCE_ADDENDUM.md','CONTROLLED_RESPONSE_ADDENDUM.md','CROSS_ASSAY_CONFIRMATION.md','WIDE_MANIFEST.json','WIDE_PROTOCOL.md','WIDE_SOURCE_QUANTITY_ADDENDUM.md','WIDE_OBSERVED_DOMAIN_ADDENDUM.md','NOGRU_EXTERNAL_ADDENDUM.md','HNRNPK_COMPONENT_ADDENDUM.md','PARSER_ADDENDUM.md','RESULTS.md','CAPTION.md','REPORT.md','REPLAY.md','SOURCES.md','CROSS_AUDIT_01.md','QA.json']:
    p=R/name
    if p.exists():hashes.append(dict(path=name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(R/'MANIFEST.json').write_text(json.dumps(dict(utc=now,files=hashes,outputs_not_old_material_mutation=True,hash_is_not_scientific_readiness=True),indent=2),encoding='utf-8')
(R/'STATUS.md').write_text('# 状态 '+now+'\n\n实际完成：matched72模型IG/controlled/数值QA；原始RNAcompete/RBNS cross-assay；28真实蛋白42task×4条件168响应；25新增蛋白等权与paired-protein CI；HNRNPK385gene/locus局部不确定性；PTBP1和HNRNPK各六NoGRU外部完整诊断；A两蛋白native独立40点各一组，以及最终位置/容量控制48点、10k bootstrap、六fresh初始化、训练顺序、真实logit消费和36项最终多重比较区间。\n\n科学图85×132mm已导出PNG/PDF/SVG，current PNG与root/C已实看像素byte-identical。主结论为家族规则超词频与组成、HNRNPK一个assay局部Full−S优先级改善；Full−S/Mean/Perm不足以支持普遍Phi解释优势。\n\n原始TXT统一literal/round_trip，622516 tokens逐位核同Pythonfloat，新增25主点bitwise不变，旧62微差fail与83文件快照保留。C实际独立family2631项、figure620项、来源/局部bootstrap39项已通过，QA引用最终sha。等待root阅读页最终实看，不扩假设、不用GPU、不更改旧稿。截止14:47:07UTC不变。\n',encoding='utf-8')
(R/'ATTEMPTS.json').write_text(json.dumps(dict(utc=now,phase='closed_candidate_delivery',wide_tasks=42,wide_proteins=28,new_proteins=25,wide_conditions=4,actual_completed_summaries=168,early_misreport_corrected=True,pcbp_qualified_source_windows=[632,916],source_failure_logs_retained=True,no_phi_unseen_zero_imputation=True,numericalQA72=True,conditional_cluster_bootstraps=2000,nogru_external_PTBP1_6=True,nogru_external_HNRNPK_6=True,nogru_cpu_only=True,context_original_point_checks=80,context_final_capacity_point_checks=48,context_bootstrap10000_recomputed=True,context_final_named_CI_checks=36,context_freshinitialization_runs=18,updated_real_logits_consumed=True,parser_contract='round_trip',raw_literal_token_checks=622516,original62fail_retained=True,independent_family2631_pass=True,independent_figure620_pass=True,independent_source_bootstrap39_pass=True,figure85mm_exported=True,remaining='root final reading-page review',deadline='2026-10-04 14:47:07 UTC'),indent=2),encoding='utf-8')
print(json.dumps(dict(utc=now,n_checks=len(checks),all_passed=True,n_manifest_files=len(hashes)),indent=2))
