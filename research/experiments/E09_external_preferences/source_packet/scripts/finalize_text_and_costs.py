"""Finalize bounded prose/replay metadata and cost rows from actual run logs."""
from pathlib import Path
import json
import pandas as pd
from literal_tables import read_csv_literal
R=Path(__file__).resolve().parents[1]
p=R/'CAPTION.md';s=p.read_text(encoding='utf-8').replace('HNRNPK RNAcompete and RBNS;', 'HNRNPK RNAcompete and the preselected 20 nM RBNS column;')
p.write_text(s,encoding='utf-8');(R/'figures'/'CAPTION.md').write_text(s,encoding='utf-8')
p=R/'CROSS_AUDIT_01.md';s=p.read_text(encoding='utf-8')
if 'Bonferroni两主比较' not in s:s+='\n最终contrast闭合：按实际`comparison`字段区分Full−control与position_capacity−NoGRU，全部36个Bonf1/2/3/15/45/90区间从已独立重建的10k draws再次核对。Bonferroni两主比较区间Full−Local+pos[−.0395639,−.0247504]，数值完全一致；首次脚本使用了消息内简写而触发assert，修为实际CSV字段后通过，未混淆两种差值。\n'
p.write_text(s,encoding='utf-8')
p=R/'REPLAY.md';s=p.read_text(encoding='utf-8')
s+='''

Final numerical contract: original author text and derived CSVs use `scripts/literal_tables.py` (`float_precision='round_trip'`), uniformly, without tolerance-based merging. Run `evaluate_rules.py`, `evaluate_wide.py`, then `replay_parser_correction.py`, `hnrnpk_stability.py`, `parser_qa.py`, `plot_external.py`, and `package_packet.py` for a frozen-array replay. The parser correction replay performs no model inference or training; archived pre-correction files and hashes are under `revisions/parser_default_20261004`.

Additional fixed HNRNPK component diagnostic: `prepare_hnr_component.py` and `evaluate_hnr_component.py` consume the saved six NoGRU count-response arrays; no new head, background or word length was chosen. All original RBNS columns and RNAcompete measurements remain separate.

Final context-control audit: `audit_position_capacity.py` is an independent CPU2 fresh initialization/implementation/epoch-order check. `cross_audit_capacity_context.py` recomputes 48 main points and all 10,000 paired gene draws from A's original profiles. `check_context_final_ci.py` explicitly distinguishes the two actual contrast names and checks 36 saved multiplicity intervals. It reads all necessary original arms and the near-capacity positional control without modifying A. The superseded pure-101-bias training plan was never run. These audit scripts require A's retained files; they do not train a model.

`results/ACTUAL_RUN_COSTS.csv` reports sums of individually recorded run elapsed seconds, not the campaign's parallel wall time. Unknown costs are not imputed. Numerical QA, source failures and the original independent 62 parser microdifference failures remain distinct from scientific outcomes.
'''
p.write_text(s,encoding='utf-8')
p=R/'SOURCES.md';s=p.read_text(encoding='utf-8')
s+='\n\nFinal text parsing contract: original score-token floats are read with round_trip, equal to Python float(token). All 622,516 selected original measured tokens were independently compared bitwise. Original authors’ scores and source archives are not edited; corrected support tables and pre-correction snapshots are separated. See PARSER_ADDENDUM.md/PARSER_QA.json.\n'
p.write_text(s,encoding='utf-8')
rows=[]
for f,phase in [('IG_AUDIT.csv','native_IG_primary72'),('CONTROLLED_AUDIT.csv','controlled_matched72'),
                 ('WIDE_EXECUTION_QA.csv','wide42task_fourconditions168'),('NOGRU_EXTERNAL_QA.csv','PTBP1_NoGRU_external6'),
                 ('HNR_COMPONENT_QA.csv','HNRNPK_NoGRU_external6')]:
    df=read_csv_literal(R/'results'/f)
    rows.append(dict(phase=phase,n_runs=len(df),summed_run_elapsed_seconds=float(df.elapsed_seconds.sum()),
                     source=f,cost_is_individual_run_sum_not_campaign_parallel_wall=True))
for f,phase in [('PARSER_REPLAY_EXECUTION.json','parser_local_frozen_replay')]:
    data=json.loads((R/'results'/f).read_text());rows.append(dict(phase=phase,n_runs=len(data),
        summed_run_elapsed_seconds=sum(x['seconds'] for x in data),source=f,cost_is_individual_run_sum_not_campaign_parallel_wall=True))
pd.DataFrame(rows).to_csv(R/'results'/'ACTUAL_RUN_COSTS.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
