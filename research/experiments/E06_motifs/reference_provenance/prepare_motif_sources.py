"""Bind the frozen S6 comparisons to exact CISBP records and stored Tomtom alignments.

No motif search, model inference, significance recalculation or Tomtom rerun.
Run using D:/Anaconda/envs/rbp/python.exe -X utf8.
"""
from pathlib import Path
import csv, hashlib, json, math
import numpy as np
import pandas as pd

WORK = Path(__file__).resolve().parent
PAPER = WORK.parent
ORIGINAL_MOTIF = Path('M:/RBP/26年数据RBP/可解释性分析/Motif')
MOTIF = WORK/'source_inputs' if (WORK/'source_inputs/outputs').exists() else ORIGINAL_MOTIF
REF = MOTIF / 'reference_data/CISBP-RNA/homo_sapiens'
S6 = PAPER / 'manuscript/source_data/S1_S14_original/Supplementary_Data_S6_motif_comparisons.csv'
S10 = PAPER / 'manuscript/source_data/S1_S14_original/Supplementary_Data_S10_motif_leaf_mapping.csv'
FIGA = Path('M:/RBP/26年数据RBP/图/figure5/figureA_data')
if (WORK/'source_inputs/Supplementary_Data_S6_motif_comparisons.csv').exists():
    S6 = WORK/'source_inputs/Supplementary_Data_S6_motif_comparisons.csv'
if (WORK/'source_inputs/Supplementary_Data_S10_motif_leaf_mapping.csv').exists():
    S10 = WORK/'source_inputs/Supplementary_Data_S10_motif_leaf_mapping.csv'
if (WORK/'source_inputs/figureA_dist_matrix.csv').exists():
    FIGA = WORK/'source_inputs'
MAIN_RBPS = ['QKI', 'PUM2', 'HNRNPA1', 'HNRNPL', 'RBFOX2', 'U2AF2']
ASSAYS = {'RNAcompete', 'HTR-SELEX', 'RNA-Bind-n-Seq', 'SEQRS'}
RULE = {
    'rbps_in_display_order': MAIN_RBPS,
    'rationale': 'Frozen named examples linking the pattern atlas, heterogeneous sequence preferences and retained PUM2/HNRNPA1 sensitivity examples; not an estimate of population validation rate.',
    'eligibility': 'Exact same RBP_Name and Motif_ID; at least one D record; measured Motif_Type in RNAcompete/HTR-SELEX/RNA-Bind-n-Seq/SEQRS; stored positive pattern 0 comparison.',
    'choice_within_rbp': 'Ascending task identifier (HepG2 before K562 where both are available), then ascending Target_ID. p/q and model scores are not sorting fields.',
    'statistics': 'Original S6 p, E and q retained verbatim; original Tomtom database/search family unchanged.',
    'main_display_window': 'Shared coordinates centered on the stored aligned reference; width max(9, reference width + 2). Matrix values and p/q unchanged.',
}

def rows(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))

def meme(path):
    out = {}
    lines = path.read_text(encoding='utf-8').splitlines()
    for i, line in enumerate(lines):
        if line.startswith('MOTIF '):
            key = line.split()[1]
            n = int(lines[i+1].split('w=')[1].split()[0])
            x = np.array([[float(v) for v in l.split()] for l in lines[i+2:i+2+n]])
            out[key] = x / x.sum(axis=1, keepdims=True)
    return out

def category(rec):
    if rec['RBP_Status'] == 'D' and rec['Motif_Type'] in ASSAYS:
        return 'direct_experiment'
    if rec['Motif_Type'] == 'JPLE':
        return 'JPLE_inferred'
    if rec['RBP_Status'] == 'I':
        return 'homology_inferred'
    if rec['RBP_Status'] == 'D':
        return 'direct_database_compiled'
    return 'unclassified'

def main():
    source = rows(S6)
    metadata = rows(REF / 'RBP_Information_all_motifs.txt', '\t')
    exact = {}
    for r in metadata:
        exact.setdefault((r['RBP_Name'], r['Motif_ID']), []).append(r)
    bound, matrices, used_meta = [], {}, []
    qcache, tcache = {}, {}
    for index, r in enumerate(source):
        key = (r['cisbp_RBP_Name'], r['Target_ID'])
        matches = exact.get(key, [])
        assert matches, ('missing exact metadata', key)
        # Keep every metadata record, including D/I ambiguity, in a separate table.
        for midx, mr in enumerate(matches):
            used_meta.append(dict(comparison_id=f'C{index+1:03}', metadata_record=midx+1, **mr))
        chosen = sorted(matches, key=lambda z: (category(z) != 'direct_experiment', z['RBP_Status'] != 'D', z['RBP_ID']))[0]
        if r['task'] not in qcache:
            qcache[r['task']] = meme(MOTIF / 'outputs' / r['task'] / 'tfmodisco_export/queries.meme')
            tp = MOTIF / 'outputs' / r['task'] / 'tomtom_cisbp/tomtom_out/tomtom.tsv'
            tcache[r['task']] = pd.read_csv(tp, sep='\t', comment='#').dropna(subset=['Query_ID', 'Target_ID'])
        tm = tcache[r['task']]
        hit = tm[(tm.Query_ID == r['Query_ID']) & (tm.Target_ID == r['Target_ID'])]
        assert len(hit) == 1, (r, len(hit))
        t = hit.iloc[0]
        for a, b in [('p_value', 'p-value'), ('q-value', 'q-value'), ('E-value', 'E-value')]:
            assert math.isclose(float(r[a]), float(t[b]), rel_tol=1e-12, abs_tol=1e-14), (r, b)
        qp = qcache[r['task']][r['Query_ID']]
        rp = pd.read_csv(REF / 'pwms' / (r['Target_ID'] + '.txt'), sep='\t')[list('ACGU')].values.astype(float)
        rp = rp / rp.sum(axis=1, keepdims=True)
        cid = f'C{index+1:03}'
        matrices[cid + '_query'] = qp
        matrices[cid + '_reference'] = rp
        out = dict(comparison_id=cid, **r)
        for field in ['RBP_Status', 'Motif_Type', 'MSource_Identifier', 'MSource_Type', 'MSource_Author', 'MSource_Year', 'PMID', 'RBP_ID']:
            out[field] = chosen[field]
        out.update(reference_class=category(chosen), metadata_n_records=len(matches),
                   reference_database='CISBP-RNA',reference_database_build='2.00',
                   metadata_all_statuses=';'.join(sorted({z['RBP_Status'] for z in matches})),
                   query_width=len(qp), reference_width=len(rp),
                   Optimal_offset=int(t['Optimal_offset']), Orientation=t['Orientation'],
                   Overlap=int(t['Overlap']), Query_consensus=t['Query_consensus'], Target_consensus=t['Target_consensus'],
                   q_lt_0_05=float(r['q-value']) < .05, q_lt_0_01=float(r['q-value']) < .01,
                   main_selected=False, main_order='', selection_rule='')
        bound.append(out)
    for order, rbp in enumerate(MAIN_RBPS):
        cand = [z for z in bound if z['cisbp_RBP_Name'] == rbp and z['reference_class'] == 'direct_experiment' and z['Query_ID'] == 'pos_pattern_0']
        assert cand, rbp
        z = sorted(cand, key=lambda z: (z['task'], z['Target_ID']))[0]
        z.update(main_selected=True, main_order=order+1, selection_rule='frozen_rbp_task_target_rule')
    for z in bound:
        rs=-z['Optimal_offset'];width=max(9,z['reference_width']+2)
        mid=rs+(z['reference_width']-1)/2
        z['main_display_alignment_min']=mid-width/2 if z['main_selected'] else ''
        z['main_display_alignment_max']=mid+width/2 if z['main_selected'] else ''
        z['main_display_nt_spacing']='shared between the query and reference' if z['main_selected'] else ''
        z['main_display_uses_aligned_window']=z['main_selected']
        z['main_query_display_trims_input']=bool(z['main_selected'] and (mid-width/2>-.5 or mid+width/2<z['query_width']-.5))
        z['main_reference_display_trims_input']=False
    for path, data in [(WORK/'motif_comparisons_all_193.csv', bound), (WORK/'CISBP_exact_metadata_records.csv', used_meta)]:
        with path.open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader(); writer.writerows(data)
    np.savez_compressed(WORK/'paired_logo_matrices.npz', **matrices)
    (WORK/'MAIN_SELECTION_RULE.json').write_text(json.dumps(RULE, indent=2), encoding='utf-8')
    # Freeze before rendering: all later rendering reads these marked rows.
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S6, S10, REF/'RBP_Information_all_motifs.txt', FIGA/'figureA_dist_matrix.csv', WORK/'MAIN_SELECTION_RULE.json']}
    d = pd.DataFrame(bound)
    summary = {
        'comparisons':len(d), 'patterns':len(d[['task','Query_ID']].drop_duplicates()),
        'tasks':d.task.nunique(), 'rbps':d.cisbp_RBP_Name.nunique(),
        'significant_comparisons':int(d.q_lt_0_05.sum()),
        'significant_patterns':len(d[d.q_lt_0_05][['task','Query_ID']].drop_duplicates()),
        'significant_tasks':d[d.q_lt_0_05].task.nunique(),
        'significant_rbps':d[d.q_lt_0_05].cisbp_RBP_Name.nunique(),
        'reference_class_counts':d.reference_class.value_counts().to_dict(),
        'reference_class_significant_counts':d[d.q_lt_0_05].reference_class.value_counts().to_dict(),
        'main_selection':d[d.main_selected].sort_values('main_order').to_dict('records'),
        'source_hashes':hashes,
        'statistical_family_note':'These q values describe each original Tomtom search; no q values recomputed after restricting to cognate/experimental references.'
    }
    assert (summary['comparisons'],summary['patterns'],summary['tasks'],summary['rbps']) == (193,80,45,31)
    assert (summary['significant_comparisons'],summary['significant_patterns'],summary['significant_tasks'],summary['significant_rbps']) == (72,35,19,13)
    (WORK/'SOURCE_AUDIT.json').write_text(json.dumps(summary, indent=2, default=int), encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['source_hashes','main_selection']}, indent=2))
    print(d[d.main_selected][['cisbp_RBP_Name','task','Target_ID','p_value','q-value','reference_class']].to_string(index=False))

if __name__ == '__main__':
    main()
