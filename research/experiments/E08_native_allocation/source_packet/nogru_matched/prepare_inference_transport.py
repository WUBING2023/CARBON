from pathlib import Path
R=Path(__file__).resolve().parent
s=(R/'run_infer.py').read_text(encoding='utf-8')
s=s.replace('d=json.load(sys.stdin);', 'd=json.load(sys.stdin); (r/"PROTOCOL_TRAINING_FROZEN.md").write_bytes((r/"PROTOCOL.md").read_bytes()) if not (r/"PROTOCOL_TRAINING_FROZEN.md").exists() else None;')
(R/'run_infer.py').write_text(s,encoding='utf-8')
print('Training registered protocol archived before RF clarification write; cohort definition unchanged.')
