"""Rerun the verified configuration without overwriting prior evidence."""
import json
from datetime import datetime
from recover_w2l import HERE,run_child

if __name__=='__main__':
    config=json.loads((HERE/'w2l_verified_config.json').read_text())
    name='verified_'+datetime.now().strftime('%Y%m%d_%H%M%S')
    run_child(dict(name=name,knobs=config['knobs'],backend='exact_fast'))
