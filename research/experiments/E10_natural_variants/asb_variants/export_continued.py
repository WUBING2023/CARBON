from pathlib import Path
import fitz
R=Path(__file__).resolve().parent;F=R/'figures';pt=72/25.4
doc=fitz.open(F/'FIG5_ADDITIONAL_a_e_85mm.pdf');page=doc[0]
span=next(s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text']=='a')
page.add_redact_annot(fitz.Rect(span['bbox']),fill=(1,1,1));page.apply_redactions(images=0,graphics=0)
from matplotlib import font_manager
fontfile=font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'))
page.insert_font(fontname='ArialBold',fontfile=fontfile)
page.insert_text(fitz.Point(*span['origin']),'d',fontsize=span['size'],fontname='ArialBold')
doc.save(F/'FIG5_ADDITIONAL_d_e_85mm.pdf');page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(F/'FIG5_ADDITIONAL_d_e_85mm.png');doc.close()
svg=(F/'FIG5_ADDITIONAL_a_e_85mm.svg').read_text('utf8').replace('>a</text>','>d</text>');(F/'FIG5_ADDITIONAL_d_e_85mm.svg').write_text(svg,'utf8')
print('Fig5 continued: original a/b/c unchanged; new natural d and cohort/conservation e; 85x91mm.')
