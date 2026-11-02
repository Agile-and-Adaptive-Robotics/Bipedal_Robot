"""Static audit of saved execution ledgers. Never imports a model or runs physics."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
BASE = ROOT / 'Neuromechanical_Models/Mujoco_SNS_models'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def table(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in r.items()})

def main():
    summary = []
    manifests = []
    for family in ('s3k', 'w2lvar', 'syn6', 'W2L704_80'):
        folder = BASE / ('W2L_704_80_baseline' if family == 'W2L704_80' else 'Gait2392_standing_families')
        source = folder / ('W2L_704_80_graph_ledger.json' if family == 'W2L704_80' else family + '_graph_ledger.json')
        d = json.loads(source.read_text(encoding='utf-8'))
        dest = OUT / family
        dest.mkdir(exist_ok=True)
        nodes, edges = d['nodes'], d['chemical_edges']
        ids = [n['id'] for n in nodes]
        errors = []
        if len(ids) != len(set(ids)):
            errors.append('duplicate node IDs')
        for e in edges:
            for field in ('source', 'target'):
                index = e[field + '_index']
                if not 0 <= index < len(nodes) or nodes[index]['name'] != e[field]:
                    errors.append({'edge': e['id'], 'bad_endpoint': field})
        table(dest / 'nodes.csv', nodes)
        table(dest / 'edges.csv', edges)
        table(dest / 'zero_pairs.csv', d.get('disabled_source_pairs', []))
        inputs = d.get('input_adapters', d.get('physical_inputs', []))
        table(dest / 'inputs.csv', inputs)
        # Enumerate each input's reachable edges once; no connectivity inferred from labels.
        adjacency = {}
        for e in edges:
            adjacency.setdefault(e['source'], []).append(e)
        paths = []
        for inp in inputs:
            start = inp['target']
            seen, queue = set(), [start]
            reached_edges = []
            while queue:
                node = queue.pop(0)
                if node in seen:
                    continue
                seen.add(node)
                for e in adjacency.get(node, []):
                    reached_edges.append(e['id'])
                    queue.append(e['target'])
            paths.append({'input_target': start, 'reachable_edge_ids': reached_edges, 'reachable_node_count': len(seen), 'interpretation': 'static reachability only; no physiological or causal validation'})
        table(dest / 'input_reachability.csv', paths)
        checks = []
        for name, c in d.get('source_checks', {}).items():
            p = ROOT / 'Code/MuJoCo_SNS/spinal' / name
            checks.append({'file': str(p.relative_to(ROOT)), 'expected': c['saved'], 'actual': digest(p) if p.exists() else None, 'matches_saved': p.exists() and digest(p) == c['saved']})
        # One byte-exact, editable saved graph checkpoint; no winner-labelled copies.
        graph = folder / ('W2L_704_80_exact_connectome.json' if family == 'W2L704_80' else family + '/exact_connectome.json')
        checkpoint = dest / 'editable_checkpoint.json'
        checkpoint.write_bytes(graph.read_bytes())
        overlay = {k: v for k, v in d.items() if k not in ('nodes', 'chemical_edges', 'disabled_source_pairs', 'compiled_edge_matrix')}
        (dest / 'runtime_context.json').write_text(json.dumps(overlay, indent=2), encoding='utf-8')
        result = {'family': family, 'ledger': str(source.relative_to(ROOT)), 'ledger_sha256': digest(source), 'census': d['census'], 'nodes_exported': len(nodes), 'edges_exported': len(edges), 'zero_pairs_exported': len(d.get('disabled_source_pairs', [])), 'endpoint_errors': errors, 'source_checks': checks, 'checkpoint_sha256': digest(checkpoint), 'scope': d.get('scope'), 'review_complete': False}
        summary.append(result)
        manifest = json.loads((folder / 'SOURCE_MANIFEST.json').read_text())
        records = manifest['files']
        if isinstance(records, dict):
            records = [{'path': k, 'sha256': v} for k, v in records.items()]
        for r in records:
            p = folder / r['path']
            manifests.append({'family': family, 'file': str(p.relative_to(ROOT)), 'exists': p.exists(), 'matches': p.exists() and digest(p) == r['sha256']})
    (OUT / 'STATIC_VALIDATION.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    (OUT / 'SOURCE_MANIFEST_CHECK.json').write_text(json.dumps(manifests, indent=2), encoding='utf-8')
    table(OUT / 'OUTPUT_MANIFEST.csv', [{'path': str(p.relative_to(OUT)), 'sha256': digest(p)} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name != 'OUTPUT_MANIFEST.csv'])
    print(json.dumps({'families': len(summary), 'endpoint_errors': sum(len(s['endpoint_errors']) for s in summary), 'manifest_mismatches': sum(not r['matches'] for r in manifests), 'source_hash_mismatches': sum(not r['matches_saved'] for s in summary for r in s['source_checks'])}))

if __name__ == '__main__':
    main()
