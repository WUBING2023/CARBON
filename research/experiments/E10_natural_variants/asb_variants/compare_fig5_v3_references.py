"""Whole original references / pre-revision / final panels, each exactly 178 mm wide."""
from pathlib import Path
import json,hashlib,shutil
import fitz
from PIL import Image
R=Path(__file__).resolve().parent;O=R/'fig5_v3';D=O/'references';D.mkdir(exist_ok=True);P=72/25.4
L=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924')
lit=L/'interpretability_strengthening_20261003/01_literature/sources'
rf=L/'supervised_expansion_20261004/candidate_portfolio_20261006/reference_faithful_20261007'
refs=[
 (lit/'prismnet_41422_2021_476_Fig7_HTML.png','PrismNet (2021), original Fig. 7','https://www.nature.com/articles/s41422-021-00476-y/figures/7','Allele identity linked to common coordinates and inputs. CARBON uses fixed reference Phi and observed BEAPR ordering; no mutant-structure or splicing-validation claim.'),
 (lit/'prismnet_41422_2021_476_Fig6_HTML.png','PrismNet (2021), original Fig. 6','https://www.nature.com/articles/s41422-021-00476-y/figures/6','Matched identity and comparative scope. CARBON retains the full cohort; no random-site/HAR control, disease enrichment or odds-ratio result is invented.'),
 (rf/'B3/references/Matreyek2018_Fig2_full.jpg','Matreyek et al. (2018), original Fig. 2','https://www.nature.com/articles/s41588-018-0122-z','Complete object matrix and aligned annotation tracks. Here, cells are frozen RBP predictions and external word preference, not abundance or pathogenicity assays.'),
 (rf/'B4/references/scrambler_Fig3.jpg','Linder et al. (2022), Scrambler Fig. 3','https://www.nature.com/articles/s42256-021-00428-6/figures/3','Literal operation identity carried through paired outputs. Here, a frozen edit/allele is matched across responses; APA function and attribution truth do not transfer.')]
v2=fitz.open(R/'fig5_v2/fig5.pdf');v3=fitz.open(O/'fig5.pdf');doc=fitz.open();metadata=[]
for i,(source,title,url,job)in enumerate(refs):
    target=D/source.name;shutil.copyfile(source,target)
    with Image.open(source)as im:iw,ih=im.size
    height=178*ih/iw
    page=doc.new_page(width=554*P,height=225*P)
    for x,title2 in[(5,'Published original: '+title),(188,'Before revision: CARBON Fig. 5 v2'),(371,'After revision: CARBON Fig. 5 v3')]:page.insert_text((x*P,8*P),title2,fontname='hebo',fontsize=10)
    page.insert_image(fitz.Rect(5*P,16*P,183*P,(16+height)*P),filename=str(target),keep_proportion=True)
    page.show_pdf_page(fitz.Rect(188*P,16*P,366*P,(16+183)*P),v2,0)
    page.show_pdf_page(fitz.Rect(371*P,16*P,549*P,(16+190)*P),v3,0)
    for x,body in[(5,job),(188,'Whole pre-revision figure. Identical width: 178 mm. All panels shown; no crop or width inflation.'),(371,'Whole final figure. Identical width: 178 mm. All 25 proteins, 36 paired tasks and 396 natural allele events retained. Negative comparative results remain visible.')]:page.insert_textbox(fitz.Rect(x*P,211*P,(x+178)*P,224*P),body,fontsize=8,fontname='helv')
    if url.startswith('https:'):page.insert_link({'kind':fitz.LINK_URI,'from':fitz.Rect(5*P,1*P,183*P,11*P),'uri':url})
    if i<2:page.get_pixmap(matrix=fitz.Matrix(120/72,120/72),alpha=False).save(O/('REFERENCE_CURRENT_V3_'+str(i+1)+'.png'))
    metadata.append({'reference':title,'original_file':str(source),'delivered_exact_original':str(target),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'original_pixel_size':[iw,ih],'comparison_width_mm':178,'comparison_height_mm':height,'source':url,'scientific_job_and_limit':job,'whole_original_not_cropped':True,'viewed_before_v3_drawing':True})
doc.save(O/'REFERENCE_CURRENT_V3_FULL_WIDTH.pdf',garbage=4,deflate=True);doc.close();v2.close();v3.close()
(O/'REFERENCE_PROVENANCE.json').write_text(json.dumps(metadata,indent=2)+'\n','utf8')
print('Four exact-original / v2 / v3 full-figure comparison pages; each panel 178 mm wide.')
