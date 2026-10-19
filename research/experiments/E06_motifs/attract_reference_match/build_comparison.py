"""Place complete reference and complete CARBON simultaneously at equal width."""
from pathlib import Path
import fitz, json, hashlib, csv
from PIL import Image

ROOT=Path(__file__).resolve().parent
P=ROOT.parent
MM=72/25.4
OUT=ROOT/'delivery';OUT.mkdir(exist_ok=True)
REF=ROOT/'references/BEAPR_Yang2019_Figure5_complete.png'
TRIAL=ROOT/'layout_trial/CARBON_BEAPR_layout_trial.pdf'
V25=ROOT/'baseline/CARBON_Fig5_v25_candidate.pdf'
URL='https://www.nature.com/articles/s41467-019-09292-w/figures/5'
def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def point(x,y):return fitz.Point(x*MM,y*MM)
def rect(x,y,w,h):return fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM)
def text(pg,x,y,s,sz=9,color=(.16,.20,.23),cn=False,bold=False):
    pg.insert_text(point(x,y),s,fontname='CN' if cn else ('hebo' if bold else 'helv'),fontsize=sz,color=color)
def lines(pg,x,y,items,sz=9,cn=False):
    for i,s in enumerate(items):text(pg,x,y+i*5,s,sz,cn=cn)

doc=fitz.open()
placements=[]
rw,rh=Image.open(REF).size
ref_h=178*rh/rw
for i,(candidate,cnname,enname) in enumerate([
    (TRIAL,'按参考结构试排','New layout trial'),
    (V25,'修改前 v25','Previous v25 candidate')]):
    pg=doc.new_page(width=420*MM,height=297*MM)
    pg.insert_font(fontname='CN',fontfile='C:/Windows/Fonts/msyh.ttc')
    text(pg,17,13,'参考原图与 CARBON：同页、同宽对照',18,cn=True)
    text(pg,17,21,'Yang et al. | Nature Communications (2019) | Figure 5 | DOI: 10.1038/s41467-019-09292-w',9.5)
    text(pg,17,29,'左：完整参考原图',12,cn=True)
    text(pg,225,29,'右：CARBON '+cnname,12,cn=True)
    text(pg,17,35,'Complete original | 178 mm wide',8,color=(.42,.47,.52))
    text(pg,225,35,enname+' | 178 mm wide',8,color=(.42,.47,.52))
    # show_pdf_page is vector-preserving; native aspect ratio, no crop.
    rd=rect(17,41,178,ref_h)
    pg.insert_image(rd,filename=str(REF),keep_proportion=True)
    with fitz.open(candidate) as src:
        pr=src[0].rect
        ch=178*pr.height/pr.width
        cd=rect(225,41,178,ch)
        pg.show_pdf_page(cd,src,0,keep_proportion=True)
    for x in [17,225]:pg.draw_line(point(x,229),point(x+178,229),color=(.80,.84,.87),width=.5)
    lines(pg,17,237,[
        '参考原图的结构：',
        '左上：群体证据；左下：位点与等位基因读出。',
        '右侧：贯穿全图、保留全部对象名称的完整比较。',
        '原图含敲低、TCGA/GTEx 与报告实验；这些结论不转移给 CARBON。'],9.5,cn=True)
    if i==0:
        notes=[
            '试排中的对应关系：',
            '左上：396 个 ASB 事件与 8 个临床注释位点。',
            '左下：ACVR1B / SLC39A13，保留相反的模型表现。',
            '右侧：全 25 个蛋白、两个固定终点；没有虚构逐蛋白分布。']
    else:
        notes=[
            '用于比较修改前后的整体结构：',
            'v25 保留 QKI 景观与原始全部数据。',
            '本页展示旧版在模块大小、右侧长列和病例联动上的差异。',
            '所有图完整显示；两侧宽度相同，未裁切或拉伸。']
    lines(pg,225,237,notes,9.5,cn=True)
    text(pg,17,266,'本次为参考与版式试排；正式正文、原候选和 Overleaf 均未替换。',10,cn=True,bold=True)
    text(pg,17,274,'Original figure: '+URL,8,color=(.23,.33,.48))
    pg.insert_link({'kind':fitz.LINK_URI,'from':rect(17,270.5,250,5),'uri':URL})
    text(pg,17,283,'Source mappings, frozen values and full cohort records accompany this PDF.',8,color=(.42,.47,.52))
    text(pg,393,283,f'{i+1} / 2',8,color=(.42,.47,.52))
    placements.append(dict(page=i+1,reference_full_image=str(REF),reference_pixel_size=[rw,rh],
                           reference_width_mm=178,reference_height_mm=ref_h,
                           CARBON_pdf=str(candidate),CARBON_width_mm=178,CARBON_height_mm=ch,
                           both_complete=True,cropped=False))
pdf=OUT/'CARBON_Fig5_reference_left_ours_right.pdf'
doc.subset_fonts()
doc.save(pdf,garbage=4,deflate=True)
doc.close()
doc=fitz.open(pdf)
for n in [0,1]:doc[n].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(OUT/f'comparison_page_{n+1}.png')
doc.close()
for ext in ['pdf','png','svg']:
    dest=OUT/f'CARBON_BEAPR_layout_trial.{ext}'
    dest.write_bytes((ROOT/f'layout_trial/CARBON_BEAPR_layout_trial.{ext}').read_bytes())
state=dict(reference={'paper':'Yang et al. Nature Communications 10, 1338 (2019)',
    'doi':'10.1038/s41467-019-09292-w','figure':'5','url':URL,
    'image':str(REF.relative_to(ROOT)),'sha256':sh(REF),'image_inspected_before_drawing':True},
    pdf=str(pdf.relative_to(ROOT)),pages=2,page_mm=[420,297],placements=placements,
    trial_native_mm=[178,151],new_experiments=False,formal_manuscript_changed=False,overleaf_changed=False)
(ROOT/'REFERENCE_SELECTION.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'pdf':str(pdf),'pages':2,'reference_size_mm':[178,ref_h],'trial_size_mm':[178,151]}))
