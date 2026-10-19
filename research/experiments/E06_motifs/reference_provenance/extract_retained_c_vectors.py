"""Retain panel-c geometry from the accepted vector PDF, with reversible axis mapping.

This is a drawing-coordinate extraction, not re-inference of model predictions.
Run in the base Python environment (PyMuPDF available).
"""
from pathlib import Path
import hashlib,json,csv
import fitz

W = Path(__file__).resolve().parent
BASE = W/'fig3_baseline.pdf'
# Axis positions read from retained vector ticks; linear transformations only.
PANELS = [
    dict(task='PUM2_K562_ENCSR661ICQ',locus='chr15:59124889–59124990 (+)',left=14.0101003647,right=199.750106812,
         ic_zero=581.326232910,ic_one=(581.326232910-540.911804199)/2,ig_zero=619.134338379,ig_one=9.243408203,
         signal_zero=666.082824707,signal_one=666.082824707-627.866821289,
         ic_range=[538,583.5],ig_range=[584,624.5], shade=[115.15599823,133.54598999],mutation=123.431396484),
    dict(task='HNRNPA1_K562_ENCSR154HRN',locus='chr1:246278962–246279063 (−)',left=209.292602539,right=393.482604980,
         ic_zero=582.083129883,ic_one=(582.083129883-541.314025879)/2,ig_zero=616.592529297,ig_one=12.38684082,
         signal_zero=666.145202637,signal_one=666.145202637-628.235351562,
         ic_range=[539,584],ig_range=[584.5,625], shade=[333.302001953,346.067993164],mutation=339.684509277),
    dict(task='HNRNPA1_HepG2_ENCSR769UEW',locus='chr8:61700784–61700885 (−)',left=405.933807373,right=588.120788574,
         ic_zero=583.061035156,ic_one=(583.061035156-543.419799805)/2,ig_zero=610.360412598,ig_one=21.431152344,
         signal_zero=666.128540039,signal_one=666.128540039-628.7113037109,
         ic_range=[540,584],ig_range=[584.5,628], shade=[436.598999023,449.226013184],mutation=442.912414551),
]
OLD_BASE = {(0.,.627,.529):'A',(.302,.733,.835):'C',(.494,.38,.282):'G',(.902,.294,.208):'U'}
SIGNALS = {(0.569,.82,.761):'eCLIP',(.902,.294,.208):'Pred. WT',(.235,.329,.533):'Pred. mut.'}

def round3(c): return tuple(round(z,3) for z in c)

def main():
    doc = fitz.open(BASE); drawings=doc[0].get_drawings()
    for cfg in PANELS:
        cfg['glyphs']={'ic':[],'ig':[]};cfg['curves']=[]
        def convert(p,track):
            n = 50 if track=='ic' else 101
            return [(.5 + n*(p.x-cfg['left'])/(cfg['right']-cfg['left'])),
                    (cfg[track+'_zero']-p.y)/cfg[track+'_one']]
        for s in drawings:
            rc=s['rect']
            if rc.x0 < cfg['left']-.15 or rc.x1 > cfg['right']+.15: continue
            if s['fill'] is not None and round3(s['fill']) in OLD_BASE:
                for tr in ['ic','ig']:
                    lo,hi=cfg[tr+'_range']
                    if rc.y0>=lo and rc.y1<=hi:
                        it=[]
                        for v in s['items']:
                            if v[0] in ['l','c']:
                                it.append([v[0]]+[convert(p,tr) for p in v[1:]])
                        if it:
                            cfg['glyphs'][tr].append(dict(base=OLD_BASE[round3(s['fill'])], items=it,closePath=s['closePath']))
            if s['color'] and round3(s['color']) in SIGNALS and len(s['items'])==100 and rc.y0>626:
                p=[s['items'][0][1]]+[v[-1] for v in s['items']]
                assert len(p)==101
                cfg['curves'].append(dict(name=SIGNALS[round3(s['color'])], points=[convert(v,'signal') for v in p]))
        assert len(cfg['curves'])==3,cfg['task']
        cfg['shade_positions']=[.5+101*(p-cfg['left'])/(cfg['right']-cfg['left']) for p in cfg['shade']]
        cfg['mutation_position']=.5+101*(cfg['mutation']-cfg['left'])/(cfg['right']-cfg['left'])
    payload=dict(source=str(BASE),source_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest(),
                 method='Linear extraction of existing vector geometry using original axis tick coordinates; no recomputation of attribution, motifs or predictions.',panels=PANELS)
    (W/'retained_c_geometry.json').write_text(json.dumps(payload),encoding='utf-8')
    with (W/'retained_c_signal_curves_from_pdf.csv').open('w',encoding='utf-8',newline='') as f:
        wr=csv.writer(f);wr.writerow(['task','curve','position_1based','normalized_intensity','source_geometry_note'])
        for c in PANELS:
            for cr in c['curves']:
                for i,(x,y) in enumerate(cr['points']):
                    assert abs(x-(i+1))<.005,(x,i)
                    wr.writerow([c['task'],cr['name'],i+1,y,'Retained PDF drawing coordinate; rounded rendering precision, not new prediction array'])
    print([(z['task'],len(z['glyphs']['ic']),len(z['glyphs']['ig']),len(z['curves'])) for z in PANELS])

if __name__=='__main__':main()
