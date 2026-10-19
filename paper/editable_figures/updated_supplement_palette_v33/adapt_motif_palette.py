"""Recolor only the retained motif-logo vector glyphs; all other raw content tokens remain."""
from pathlib import Path
import fitz,hashlib,json,re
from pypdf import PdfReader
from pypdf.generic import ContentStream
import numpy as np
O=Path(__file__).resolve().parent
I=O/'inputs'
palette=json.loads((O/'PALETTE_V3.json').read_text())['nucleotides']
old={'A':'#00A087','C':'#4DBBD5','G':'#7E6148','U':'#E64B35'}
rgb=lambda s:tuple(int(s[i:i+2],16) for i in (1,3,5))
mapping={rgb(v):rgb(palette[k]) for k,v in old.items()}
reports=[]
num=rb'[+-]?(?:\d+(?:\.\d*)?|\.\d+)'
pat=re.compile(rb'(?<![\w.])('+num+rb')\s+('+num+rb')\s+('+num+rb')\s+(rg|RG)\b')
for stem,logo_height in [('figS9_conservation_complete',31.1811023622),('figS12_previous',25.5118110236)]:
    original=I/(stem+'.pdf')
    assert original.is_file(), original
    rd=PdfReader(original);cs=ContentStream(rd.pages[0].get_contents(),rd)
    replacements={}; ordinal=0; stack=[];clip=None;lastre=None
    for args,op in cs.operations:
        if op==b'q':stack.append((clip,lastre))
        elif op==b'Q':clip,lastre=stack.pop()
        elif op==b're':lastre=[float(v) for v in args]
        elif op in (b'W',b'W*'):clip=lastre
        elif op in (b'rg',b'RG'):
            color=tuple(round(float(v)*255) for v in args)
            if clip and abs(clip[3]-logo_height)<1e-7 and clip[0]>200 and color in mapping:
                replacements[ordinal]=mapping[color]
            ordinal+=1
    pdf=fitz.open(original);streams=pdf[0].get_contents();assert len(streams)==1
    raw=pdf.xref_stream(streams[0]);matches=list(pat.finditer(raw));assert len(matches)==ordinal
    out=[];end=0; changes=[]
    for n,match in enumerate(matches):
        if n not in replacements:continue
        out.append(raw[end:match.start()]);new=replacements[n]
        val=(' '.join(f'{v/255:.10f}' for v in new)+' '+match.group(4).decode()).encode()
        out.append(val);end=match.end();changes.append({'rgb_ordinal':n,'old_token':match.group().decode(),'new_token':val.decode()})
    out.append(raw[end:]);newraw=b''.join(out);assert changes
    # Every untouched raw byte is retained; changes are confined to selected RGB operands.
    pdf.update_stream(streams[0],newraw)
    dest=O/(stem+'.pdf');pdf.save(dest,garbage=4,deflate=True);pdf.close()
    before=fitz.open(original);after=fitz.open(dest)
    assert before[0].get_text()==after[0].get_text()
    a=before[0].get_pixmap(dpi=180);b=after[0].get_pixmap(dpi=180)
    aa=np.frombuffer(a.samples,np.uint8).reshape(a.height,a.width,3);bb=np.frombuffer(b.samples,np.uint8).reshape(b.height,b.width,3)
    diff=np.any(aa!=bb,axis=2); yy,xx=np.where(diff)
    assert len(xx)>0
    bbox=[float(xx.min()*72/180),float(yy.min()*72/180),float((xx.max()+1)*72/180),float((yy.max()+1)*72/180)]
    # Only top-right logo rows may change; the adjacent CWM, heatmaps and tracks are untouched.
    assert bbox[0]>200 and bbox[3]<100,bbox
    b.save(O/(stem+'_after.png'))
    reports.append({'figure':stem,'input_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
        'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'changed_RGB_commands':changes,
        'all_text_unchanged':True,'noncolor_content_bytes_retained':True,'pixel_change_bbox_pt':bbox})
    before.close();after.close()


(O/'ADAPTER_CHECK.json').write_text(json.dumps({'repairs':reports},indent=2),encoding='utf8')
print('S9/S12 motif RGB-only vector adapters reproduced.')
