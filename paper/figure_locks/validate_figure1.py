"""Check the exact user-protected Figure 1 before compiling or packaging."""
from pathlib import Path
import hashlib,json
M=Path(__file__).resolve().parents[1]
lock=json.loads((M/'FIGURE_LOCKS.json').read_text('utf8'))['Figure1']
for path_key,hash_key in [('pdf','pdf_sha256'),('source','source_sha256')]:
    path=M/lock[path_key]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual==lock[hash_key],f'User-protected Figure 1 changed: {lock[path_key]}'
assert r'\includegraphics[width=85mm]{fig1.pdf}' in (M/'CARBON_main.tex').read_text('utf8')
print(json.dumps({'figure1_locked':True,'accepted_version':lock['accepted_version'],
                  'pdf_sha256':lock['pdf_sha256'],'placement_width_mm':85}))
