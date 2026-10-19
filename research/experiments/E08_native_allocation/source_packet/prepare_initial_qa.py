from pathlib import Path
O=Path(__file__).resolve().parent/'nogru_position_capacity'
s=(O/'run_remote.py').read_text().replace("'remote_train.py'","'remote_initial_qa.py'").replace('r/"remote_train.py"','r/"remote_initial_qa.py"').replace('REMOTE_LOG.jsonl','INITIAL_QA_LOG.json').replace('TRANSPORT_RUN.json','TRANSPORT_INITIAL_QA.json')
(O/'run_initial_qa.py').write_text(s,encoding='utf-8')
