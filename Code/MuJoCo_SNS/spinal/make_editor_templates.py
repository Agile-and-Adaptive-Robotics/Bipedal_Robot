#!/usr/bin/env python
"""Build connectome_block_editor template library (2026-09-23).

Sources (all read-only):
  - Neuromechanical_Models\\Biped_2xCPG_wSubs\\connectivity_audit.txt (261-edge
    wiring inventory)                -> template 'biped2xcpg'
  - %TEMP%\\tpl_w2l.json (subagent mine of Walker_2_Layer_CPG.aproj), fallback
    spinal\\w2l_equivalent_draft.json -> 'w2l'
  - BilateralRG documented transform of w2l (tools/SESSION_NOTES_20260916.md)
                                      -> 'bilateralrg'
  - %TEMP%\\tpl_li.json (subagent mine of Li Model walk tester rearranged.aproj)
                                      -> 'li'
  - spinal\\best_walk_params_v9/v10.json + curriculum_stage1/2/3.json +
    reports_20260923\\s3k_trial34_full_params.json -> walker_v9/v10/s1/s2/s3/s3k
  - spinal\\synergy_basis.npz (rank-6 NMF per leg) -> synergy6
  - CONNECTOME.md rule table          -> 'rules'
  - spinal\\replication\\{shevtsova,shinohara,rybak}_rules.json -> converted

Output: spinal\\connectome_templates.json (sidecar loaded by the editor).
Types are validated against the editor catalog. Re-run after any source
changes; the editor reads the file fresh on every page load.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
NM = r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
TMP = os.environ.get('TEMP', r'C:\Users\Ben Bolen\AppData\Local\Temp')

TYPES = set("""SN-Ia SN-II SN-Ib SN-heel SN-toe PORT-load IN-V0D IN-V0V IN-V1
IN-V2a IN-V2b IN-V3 IN-C HC-RG-E HC-RG-F IN-InE IN-InF HC-PF-E HC-PF-F IN-PF
IN-IaIN IN-IbIN IN-IIe IN-IIi IN-Ib+ IN-LBIN IN-KINH RC MN MUSCLE""".split())

REPLICATION_TYPEMAP = {
    'DRIVE': 'PORT-load', 'IN-inh': 'IN-C', 'CIN-exh': 'IN-V3',
    'CIN-inh': 'IN-V0D', 'V2a': 'IN-V2a', 'IaIN': 'IN-IaIN', 'RC': 'RC',
    'MN': 'MN', 'HC-RG-E': 'HC-RG-E', 'HC-RG-F': 'HC-RG-F',
    'HC-PF-E': 'HC-PF-E', 'HC-PF-F': 'HC-PF-F',
    'AFF-in': 'SN-Ia',   # refined per-label below
}

def strip_paren(l):
    return re.sub(r'\s*\([^)]*\)\s*$', '', l).strip()

def resolve_alias(spec):
    """Nodes may carry ' (...) ' suffixes the edges don't; remap edge
    endpoints to the label's prefix form and strip the suffix."""
    N = spec['nodes']
    clean2full = {}
    for n in N:
        c = strip_paren(n['label'])
        clean2full.setdefault(c, n['label'])
        n['label'] = c
    # de-duplicate labels after stripping (keep first)
    seen, drop = set(), []
    for n in N:
        if n['label'] in seen:
            drop.append(n)
        else:
            seen.add(n['label'])
    for n in drop:
        N.remove(n)
    E = spec.get('edges', spec.get('synapses', []))
    for e in E:
        for k in ('from', 'to'):
            v = strip_paren(e[k])
            if v in seen:
                e[k] = v
    spec['edges'] = E
    return spec

def aff_type(label):
    l = label.lower()
    if 'ii' in l.replace('ii_', ' '): pass
    if re.search(r'\bii\b|_ii_|ii$|-ii', l): return 'SN-II'
    if re.search(r'\bib\b|_ib_|ib$|-ib', l): return 'SN-Ib'
    return 'SN-Ia'

def convert_replication(path):
    d = json.load(open(path, encoding='utf-8'))
    nodes, edges = [], []
    for n in d.get('nodes', []):
        t = REPLICATION_TYPEMAP.get(n['type'])
        if t is None:
            t = 'IN-C'
        if n['type'] == 'AFF-in':
            t = aff_type(n['label'])
        nodes.append({'type': t, 'label': n['label'],
                      'x': n.get('x', 100), 'y': n.get('y', 100)})
    for e in d.get('synapses', d.get('edges', [])):
        edges.append({'from': strip_paren(e['from']), 'to': strip_paren(e['to']),
                      'sign': e['sign'],
                      'gain': e.get('gain', 0.5),
                      'tag': (e.get('tag') or 'lit')[:24]})
    spec = {'nodes': nodes, 'edges': edges}
    return resolve_alias(spec)


# ---------------------------------------------------------------- biped2xcpg
def build_biped2xcpg():
    """Parse connectivity_audit.txt. Signs: the audit carries gains only, so
    polarity is assigned by wiring class (Deng/AnimatLab conventions of this
    project family); classes not derivable are tagged sign_assumed."""
    txt = open(os.path.join(NM, 'Biped_2xCPG_wSubs', 'connectivity_audit.txt'),
               encoding='utf-8', errors='replace').read()
    rows = re.findall(r'^\s*\d+\.\s+(.*?)\s+=>\s+(.*?)\s+G=([\d.]+)\s*$',
                      txt, re.M)
    def clean(s):
        return s.replace('OffPage:', '').strip()
    pairs, labels = [], set()
    for a, b, g in rows:
        a, b = clean(a), clean(b)
        pairs.append((a, b, float(g)))
        labels.add(a); labels.add(b)

    def base(l):     # strip trailing spaces variants used in the audit
        return re.sub(r'\s+', ' ', l).strip()

    def cls(l):
        s = base(l)
        if re.search(r'\bRG\b', s) and 'IN' not in s:
            return 'rg'
        if 'IN RG' in s: return 'rg_in'
        if re.search(r'\bIa\b', s) and re.search(r' (E|F|2)$', s): return 'ia_in'
        if re.search(r'\b(RE|R [EF])\b', s) or re.search(r'\bR [EF]$', s): return 'rc'
        if re.search(r'\bIa\b', s): return 'sn_ia'
        if re.search(r'\bII\b', s): return 'sn_ii'
        if re.search(r'\bIb\b', s): return 'sn_ib'
        if 'PF' in s and 'IN' in s: return 'pf_in'
        if 'PF' in s: return 'pf'
        if 'MN' in s: return 'mn'
        return 'other'

    def side(l):
        return 'L' if re.match(r'^L[H_ ]', base(l)) else 'R'

    def ntype(l):
        s, c = base(l), cls(l)
        if c == 'rg':       return 'HC-RG-E' if re.search(r'(ext|E)\b', s.split()[-1]) or 'ext' in s.lower() else 'HC-RG-F'
        if c == 'rg_in':    return 'IN-InE' if 'E' in s.split()[-1] or 'ext' in s.lower() else 'IN-InF'
        if c == 'ia_in':    return 'IN-IaIN'
        if c == 'rc':       return 'RC'
        if c == 'sn_ia':    return 'SN-Ia'
        if c == 'sn_ii':    return 'SN-II'
        if c == 'sn_ib':    return 'SN-Ib'
        if c == 'pf_in':    return 'IN-PF'
        if c == 'pf':       return 'HC-PF-E' if re.search(r'(E|Ext|ext)', s.split()[-1]) else 'HC-PF-F'
        if c == 'mn':       return 'MN'
        return 'IN-C'

    # sign by (class_from, class_to)
    def sign_of(a, b):
        ca, cb = cls(a), cls(b)
        INH_FROM = {'rg_in', 'ia_in', 'rc'}
        if ca in INH_FROM: return 'inh', ''
        if ca == 'sn_ia' and cb == 'ia_in': return 'exc', ''
        if ca in ('sn_ia', 'sn_ii', 'sn_ib'): return 'exc', ''
        if ca == 'pf' and cb == 'mn': return 'exc', ''
        if ca == 'pf_in': return 'inh', ''
        if ca == 'mn' and cb == 'rc': return 'exc', ''
        if ca == 'rg' and cb in ('rg', 'rg_in', 'pf', 'pf_in'): return 'exc', 'sign_assumed' if cb == 'rg' else ''
        if cb == 'mn' and ca in ('ia_in', 'rc', 'pf_in'): return 'inh', ''
        if ca == 'rg' and cb == 'mn': return 'exc', ''
        return 'exc', 'sign_assumed'

    nodes, placed = [], {}
    order = sorted(labels, key=lambda l: (side(l), cls(l), base(l)))
    BAND = {'rg': 90, 'rg_in': 90, 'pf': 210, 'pf_in': 210, 'ia_in': 330,
            'mn': 450, 'rc': 570, 'sn_ia': 690, 'sn_ii': 750, 'sn_ib': 810}
    col = {}
    for l in order:
        c = cls(l)
        key = (side(l), c)
        col[key] = col.get(key, 0) + 1
        x = (60 if side(l) == 'L' else 620) + (col[key] - 1) * 78
        y = BAND.get(c, 640) + (side(l) == 'R') * 0
        # R side gets its own full column block to avoid label pileups
        if side(l) == 'R':
            x = 620 + (col[key] - 1) * 78
        nd = {'type': ntype(l), 'label': base(l), 'x': x,
              'y': BAND.get(c, 640)}
        nodes.append(nd); placed[base(l)] = nd
    edges = []
    for a, b, g in pairs:
        s, note = sign_of(a, b)
        edges.append({'from': base(a), 'to': base(b), 'sign': s, 'gain': g,
                      'tag': note or 'aproj'})
    # implied muscles per MN group (audit is neuron-only)
    mus = []
    for nd in nodes:
        if nd['type'] != 'MN': continue
        ml = nd['label'].replace('MN', 'muscle')
        mus.append({'type': 'MUSCLE', 'label': ml, 'x': nd['x'] + 40,
                    'y': nd['y'] + 90})
        edges.append({'from': nd['label'], 'to': ml, 'sign': 'exc',
                      'gain': 1.0, 'tag': 'implied_mn_to_muscle'})
    nodes.extend(mus)
    return {'nodes': nodes, 'edges': edges}


# ---------------------------------------------------------------- bilateralrg
def _side(l):
    if re.match(r'^L[ _]', l): return 'L'
    if re.match(r'^R[ _]', l): return 'R'
    m = re.search(r'[_ ]([LR])$', l)
    return m.group(1) if m else None

def build_bilateralrg(w2l):
    """Real W2L base + the documented BilateralRG build chain
    (tools/SESSION_NOTES_20260916.md):
      build_rg.pl  -> R RG half-center ipsilateral to R Hip/Knee PF,
                      the 4 crossed L-RG->R-PF connexions REMOVED,
                      Stimulus_2 (10 nA into R RG flx);
      build_comm.pl + patch_comm_types -> Shinohara c1/V3 commissurals;
      build_aff.pl -> flexor Ia/II + extensor Ib excite their HCs (0.01);
      build_contact.pl -> heel/toe contact -> ipsilateral extensor layers."""
    spec = json.loads(json.dumps(w2l))
    N, E = spec['nodes'], spec['edges']
    by = {n['label']: n for n in N}
    def add(t, l, x, y):
        if l not in by:
            nd = {'type': t, 'label': l, 'x': x, 'y': y}
            N.append(nd); by[l] = nd
        return l
    def e(f, t, s, g, tag):
        E.append({'from': f, 'to': t, 'sign': s, 'gain': g, 'tag': tag})

    # --- step 1: bilateral RG (remove crossed L-RG -> R-PF, add R half-center)
    crossed = [i for i, x in enumerate(E)
               if by.get(x['from'], {}).get('type') in ('HC-RG-E', 'HC-RG-F')
               and by.get(x['to'], {}).get('type') in ('HC-PF-E', 'HC-PF-F')
               and _side(x['from']) == 'L' and _side(x['to']) == 'R']
    for i in sorted(crossed, reverse=True):
        E.pop(i)
    l_rge = next((n['label'] for n in N if n['type'] == 'HC-RG-E'), None)
    l_rgf = next((n['label'] for n in N if n['type'] == 'HC-RG-F'), None)
    if l_rge and _side(l_rge) == 'L':
        y0 = by[l_rge]['y']
        r_rge = add('HC-RG-E', 'R RG ext', 1110, y0 + 150)
        r_rgf = add('HC-RG-F', 'R RG flx', by[l_rgf]['x'] + 1000,
                    by[l_rgf]['y'] + 150)
        add('IN-InE', 'R RG ext IN', by[l_rge]['x'] + 1000, y0 + 90)
        add('IN-InF', 'R RG flx IN', by[l_rgf]['x'] + 1000,
            by[l_rgf]['y'] + 90)
        add('PORT-load', 'Stimulus_2 (10nA)', 1110, y0 - 60)
        e('Stimulus_2 (10nA)', r_rgf, 'exc', 10.0, 'kickoff_antiphase')
        e(r_rge, 'R RG ext IN', 'exc', 0.5, 'rg_laminate')
        e('R RG ext IN', r_rgf, 'inh', 0.5, 'rg_laminate')
        e(r_rgf, 'R RG flx IN', 'exc', 0.5, 'rg_laminate')
        e('R RG flx IN', r_rge, 'inh', 0.5, 'rg_laminate')
        # ipsilateral drives mirror how L RG drives L PF
        for x in list(E):
            f = by.get(x['from'], {})
            t = by.get(x['to'], {})
            if f.get('type') in ('HC-RG-E', 'HC-RG-F') and \
               t.get('type') in ('HC-PF-E', 'HC-PF-F') and \
               _side(x['from']) == 'L' and _side(x['to']) == 'L':
                nt = next((n['label'] for n in N if n['type'] == t['type'] and
                           _side(n['label']) == 'R' and
                           n['label'].split()[1] == x['to'].split()[1]), None)
                nf = r_rge if f['type'] == 'HC-RG-E' else r_rgf
                if nt:
                    e(nf, nt, x['sign'], x['gain'], x['tag'])

    # --- step 2: commissural c1 / V3 (Shinohara 2025)
    for S in ('L', 'R'):
        O = 'R' if S == 'L' else 'L'
        rgf_s = next((n['label'] for n in N if n['type'] == 'HC-RG-F' and
                      _side(n['label']) == S), None)
        rgf_o = next((n['label'] for n in N if n['type'] == 'HC-RG-F' and
                      _side(n['label']) == O), None)
        rge_s = next((n['label'] for n in N if n['type'] == 'HC-RG-E' and
                      _side(n['label']) == S), None)
        rge_o = next((n['label'] for n in N if n['type'] == 'HC-RG-E' and
                      _side(n['label']) == O), None)
        ine_o = next((n['label'] for n in N if n['type'] == 'IN-InE' and
                      _side(n['label']) == O), None)
        c1 = add('IN-C', 'c1_' + S, 1180, 60 if S == 'L' else 140)
        v3 = add('IN-V3', 'V3_' + S, 1180, 240 if S == 'L' else 320)
        if rgf_s and rgf_o:
            e(rgf_s, c1, 'exc', 1.0, 'comm_c1')
            e(c1, rgf_o, 'inh', 2.749, 'c1_SynAmp2.749')
        if rge_s and rge_o:
            e(rge_s, v3, 'exc', 1.0, 'comm_v3')
            e(v3, rge_o, 'exc', 0.1, 'v3_SynAmp0.1_weak')
        if v3 and ine_o:
            e(v3, ine_o, 'exc', 0.1, 'v3_to_contra_InE')

    # --- step 3: afferent -> HC excitation ("Afferent HC Excite" SynAmp 0.01)
    for aff in [n for n in N if n['type'] in ('SN-Ia', 'SN-II', 'SN-Ib')]:
        if aff['type'] == 'SN-Ib':
            kind, pat = 'E', r'(ext|Ext|stance|Stance)'
        else:
            kind, pat = 'F', r'(flx|Flx|flex|swing|Swing)'
        ismatch = re.search(pat, aff['label'])
        s = _side(aff['label'])
        if not ismatch or not s:
            continue
        for tgt in N:
            if tgt['type'] == 'HC-PF-' + kind and _side(tgt['label']) == s and \
               re.search(pat if kind == 'F' else pat, tgt['label']):
                e(aff['label'], tgt['label'], 'exc', 0.01, 'aff_HC_excite')

    # --- step 4: contact -> ipsilateral extensor layers (Gain C = 20)
    for S in ('L', 'R'):
        heel = add('SN-heel', S + ' heel contact', 40 if S == 'L' else 1140, 940)
        toe = add('SN-toe', S + ' toe contact', 140 if S == 'L' else 1240, 940)
        for tgt in N:
            ts = _side(tgt['label'])
            if ts != S:
                continue
            if tgt['type'] in ('HC-RG-E', 'HC-PF-E') and \
               re.search(r'(ext|Ext|stance)', tgt['label']):
                e(heel, tgt['label'], 'exc', 20.0, 'contact_C20')
                e(toe, tgt['label'], 'exc', 20.0, 'contact_C20')
            elif tgt['type'] == 'MN' and re.search(r'(ext|Ext)', tgt['label']):
                e(heel, tgt['label'], 'exc', 20.0, 'contact_C20')
                e(toe, tgt['label'], 'exc', 20.0, 'contact_C20')
    return spec


# ---------------------------------------------------------------- walker v-series
POOLS = ['hip_ext', 'hip_flex', 'knee_ext', 'knee_flex', 'ankle_pf', 'ankle_df']

def build_walker(name, params, pf_gain=None, note=''):
    p = params
    N, E = [], []
    def nd(t, l, x, y):
        N.append({'type': t, 'label': l, 'x': x, 'y': y}); return l
    def e(f, t, s, g, tag):
        E.append({'from': f, 'to': t, 'sign': s,
                  'gain': (round(g, 3) if isinstance(g, (int, float)) else g),
                  'tag': tag})
    def g(k, d=None):
        return p.get(k, d)
    X0 = 60
    nd('PORT-load', 'DRIVE', X0, 120)
    nd('PORT-load', 'POSTURE', X0, 320)
    if pf_gain is not None:
        nd('PORT-load', 'pf_gain=%.3f' % pf_gain, X0, 200)
    for i, S in enumerate(('R', 'L')):
        X = X0 + 240 + i * 950
        rg_e = nd('HC-RG-E', 'RG_E_' + S, X, 120)
        rg_f = nd('HC-RG-F', 'RG_F_' + S, X, 260)
        ine = nd('IN-InE', 'InE_' + S, X - 120, 120)
        inf = nd('IN-InF', 'InF_' + S, X - 120, 260)
        e('DRIVE', rg_e, 'exc', g('desc_e', 1.0), 'desc_e')
        e('DRIVE', rg_f, 'exc', g('desc_f', 1.0), 'desc_f')
        e(rg_e, ine, 'exc', 1.0, 'rg_laminate')
        e(inf, rg_e, 'inh', g('rg_adapt', 1.0) or 1.0, 'rg_adapt')
        e(rg_f, inf, 'exc', 1.0, 'rg_laminate')
        e(ine, rg_f, 'inh', g('rg_adapt', 1.0) or 1.0, 'rg_adapt')
        # commissural c1 / V3
        c1 = nd('IN-C', 'c1_' + S, X + 110, 300)
        v3 = nd('IN-V3', 'V3_' + S, X + 110, 80)
        OS = 'L' if S == 'R' else 'R'
        XO = X0 + 240 + (1 - i) * 950
        e(rg_f, c1, 'exc', 1.0, 'comm_c1')
        e(c1, 'RG_F_' + OS, 'inh', g('c1_gain', 0.0), 'c1_gain')
        e(rg_e, v3, 'exc', 1.0, 'comm_v3')
        e(v3, 'RG_E_' + OS, 'exc', g('v3_gain', 0.0), 'v3_gain')
        e(v3, 'InE_' + OS, 'exc', g('v3_gain', 0.0), 'v3_to_contraInE')
        # PF groups
        pe1 = nd('HC-PF-E', 'PF_E1_' + S, X + 250, 90)
        pe2 = nd('HC-PF-E', 'PF_E2_' + S, X + 250, 160)
        pf1 = nd('HC-PF-F', 'PF_F1_' + S, X + 250, 230)
        pf2 = nd('HC-PF-F', 'PF_F2_' + S, X + 250, 300)
        e(rg_e, pe1, 'exc', g('rg_to_pf', 1.0), 'rg_to_pf')
        e(rg_e, pe2, 'exc', g('e2_pf', g('rg_to_pf', 1.0)), 'e2_pf')
        e(rg_f, pf1, 'exc', g('rg_to_pf', 1.0), 'rg_to_pf')
        e(rg_f, pf2, 'exc', g('rg_to_pf', 1.0), 'rg_to_pf')
        # MN pools + muscles
        mn = {}
        for j, pool in enumerate(POOLS):
            mn[pool] = nd('MN', 'MN_' + pool + '_' + S, X + 430, 60 + j * 58)
            nd('MUSCLE', pool.replace('_', ' ') + '_' + S,
               X + 600, 60 + j * 58)
            e(mn[pool], pool.replace('_', ' ') + '_' + S, 'exc', 1.0,
              'mn_to_muscle')
        pg = pf_gain if pf_gain is not None else 1.0
        for pool in ('hip_ext', 'knee_ext', 'ankle_pf'):
            e(pe1, mn[pool], 'exc', pg, 'pf_gain')
            e(pe2, mn[pool], 'exc', pg, 'pf_gain')
        for pool in ('hip_flex', 'knee_flex'):
            e(pf1, mn[pool], 'exc', pg, 'pf_gain')
            e(pf2, mn[pool], 'exc', pg, 'pf_gain')
        if g('f1_df') is not None:
            e(pf1, mn['ankle_df'], 'exc', g('f1_df'), 'f1_df')
        if g('f1_kf') is not None:
            e(pf1, mn['knee_flex'], 'exc', g('f1_kf'), 'f1_kf')
        if g('post_kneext') is not None:
            e('POSTURE', mn['knee_ext'], 'exc', g('post_kneext'), 'post_kneext')
        if g('post_hipext') is not None:
            e('POSTURE', mn['hip_ext'], 'exc', g('post_hipext'), 'post_hipext')
        # trim knob annotation on ankle PF posture bias
        if g('ankle_post_walk_trim') is not None:
            e('POSTURE', mn['ankle_pf'], 'exc', g('ankle_post_walk_trim'),
              'ankle_post_walk_trim')
        # KINH swing suppression
        kinh = nd('IN-KINH', 'KINH_' + S, X + 340, 380)
        e(pf1, kinh, 'exc', 1.0, 'f1_gate')
        if g('f1_kneext_inh'):
            e(kinh, mn['knee_ext'], 'inh', g('f1_kneext_inh'), 'f1_kneext_inh')
        if g('f1_anklepf_inh'):
            e(kinh, mn['ankle_pf'], 'inh', g('f1_anklepf_inh'), 'f1_anklepf_inh')
        if g('contra_kinh') is not None:
            e(pf1, 'KINH_' + OS, 'exc', g('contra_kinh'), 'contra_kinh')
        # phase-reset INs (v9/v10 era)
        if g('phase_reset_e') is not None:
            pre = nd('IN-C', 'PRESET_E_' + S, X - 230, 60)
            prf = nd('IN-C', 'PRESET_F_' + S, X - 230, 320)
            sig_e = nd('SN-Ia', 'hip_ext_sig_' + S, X0 - 0 + 10, 40)
            sig_f = nd('SN-Ia', 'hip_flex_sig_' + S, X0 + 10, 360)
            e(sig_e, pre, 'exc', 1.0, 'hip_ext_sig')
            e(sig_f, prf, 'exc', 1.0, 'hip_flex_sig')
            e(pre, rg_e, 'exc', g('phase_reset_e'), 'phase_reset_e')
            e(prf, rg_f, 'exc', g('phase_reset_f'), 'phase_reset_f')
        # afferents (per joint, 3 classes); shared IaIN/LBIN created once
        iain = None
        if g('ia_in') is not None:
            iain = nd('IN-IaIN', 'IaIN_' + S, X + 200, 470)
            e(iain, mn['hip_ext'], 'inh', g('ia_in'), 'ia_in')
            e(iain, mn['hip_flex'], 'inh', g('ia_in'), 'ia_in')
        lbin = None
        if g('ib_rge') is not None:
            lbin = nd('IN-LBIN', 'LBIN_' + S, X + 90, 540)
            e(lbin, rg_e, 'exc', g('ib_rge'), 'ib_rge')
        affs = []
        for j, joint in enumerate(('hip', 'knee', 'ankle')):
            for k, t in (('Ia', 'SN-Ia'), ('II', 'SN-II'), ('Ib', 'SN-Ib')):
                affs.append(nd(t, k + '_' + joint + '_' + S,
                               X - 170, 470 + j * 60 + (k == 'II') * 18 +
                               (k == 'Ib') * 36))
        for a in affs:
            pref = a.rsplit('_', 2)[0]; joint = a.rsplit('_', 2)[1]
            own = {'hip': ('hip_ext', 'hip_flex'), 'knee': ('knee_ext', 'knee_flex'),
                   'ankle': ('ankle_pf', 'ankle_df')}[joint]
            ext_mn, flx_mn = mn[own[0]], mn[own[1]]
            if pref == 'Ia':
                e(a, ext_mn, 'exc', 1.0, 'ia_mono')
                e(a, flx_mn, 'exc', 1.0, 'ia_mono')
                if iain:
                    e(a, iain, 'exc', 1.0, 'ia_to_iaIN')
                else:
                    e(a, flx_mn, 'inh', 0.5, 'ia_recip_direct_legacy')
                    e(a, ext_mn, 'inh', 0.5, 'ia_recip_direct_legacy')
            elif pref == 'Ib':
                if lbin:
                    e(a, lbin, 'exc', 1.0, 'ib_to_LBIN')
                else:
                    e(a, ext_mn, 'inh', 0.5, 'ib_auto_direct_legacy')
                    e(a, flx_mn, 'inh', 0.5, 'ib_auto_direct_legacy')
            else:
                if g('ii_f_central') is not None:
                    e(a, rg_f, 'exc', g('ii_f_central'), 'ii_f_central')
                    e(a, rg_e, 'exc', g('ii_e_central') or 0.0, 'ii_e_central')
                else:
                    e(a, ext_mn, 'exc', 0.3, 'ii_aff_legacy')
        # central Ia flexor knob
        if g('ia_f_central') is not None:
            e('Ia_hip_' + S, rg_f, 'exc', g('ia_f_central'), 'ia_f_central')
        # IBEXC stance-gated reversal
        if g('v3_to_ibexc') is not None:
            ibx = nd('IN-Ib+', 'IBEXC_' + S, X + 90, 600)
            e('Ib_ankle_' + S, ibx, 'exc', 1.0, 'ib_to_IBEXC')
            e(rg_e, ibx, 'exc', 1.0, 'stance_gate')
            e(ibx, mn['ankle_pf'], 'exc', 0.5, 'ib_reversal')
            e(v3, ibx, 'exc', g('v3_to_ibexc'), 'v3_to_ibexc')
        # heel/toe
        if g('heel_rge') is not None:
            heel = nd('SN-heel', 'heel_' + S, X - 260, 560)
            toe = nd('SN-toe', 'toe_' + S, X - 260, 620)
            e(heel, rg_e, 'exc', g('heel_rge'), 'heel_rge')
            e(toe, rg_e, 'exc', g('toe_rge'), 'toe_rge')
            if g('contact_onset') is not None:
                e(heel, mn['ankle_pf'], 'exc', g('contact_onset'),
                  'contact_onset')
        # Renshaw
        rc = nd('RC', 'RC_' + S, X + 430, 420)
        for pool in ('knee_ext', 'hip_ext'):
            e(mn[pool], rc, 'exc', 1.0, 'mn_to_rc')
        e(rc, mn['knee_ext'], 'inh', g('renshaw', 0.5), 'renshaw')
        e(rc, mn['hip_ext'], 'inh', g('renshaw', 0.5), 'renshaw')
    # shared trunk + phase machine
    nd('MN', 'MN_trunk', X0, 640)
    nd('MUSCLE', 'trunk', X0 + 170, 640)
    e('POSTURE', 'MN_trunk', 'exc', 0.3, 'trunk_IMU_pd')
    e('MN_trunk', 'trunk', 'exc', 1.0, 'mn_to_muscle')
    if g('pm_gain') is not None:
        pm = nd('IN-C', 'PM T=%.2f ws=%.2f' % (g('pm_T', 0), g('pm_ws', 0)),
                X0 + 430, 640)
        e(pm, 'RG_E_R', 'exc', g('pm_gain'), 'pm_gain')
        e(pm, 'RG_E_L', 'exc', g('pm_gain'), 'pm_gain')
    return {'nodes': N, 'edges': E, '_note': note}


# ---------------------------------------------------------------- synergy6
def pool_of(m):
    m = m.lower()
    def has(*ws): return any(w in m for w in ws)
    if has('glut_max', 'glut_med', 'glut_min'): return 'hip_abd' if has('glut_med', 'glut_min') else 'hip_ext'
    if has('add_mag', 'add_long', 'add_brev', 'adductor', 'grac'): return 'hip_add'
    if has('psoas', 'iliac', 'tfl'): return 'hip_flex'
    if has('semimem', 'semiten', 'bifemsh', 'bflh'): return 'knee_flex'
    if has('rect_fem'): return 'knee_ext_hip_flex'
    if has('vas_int', 'vas_lat', 'vas_med', 'vast'): return 'knee_ext'
    if has('gastroc', 'gas_med', 'gas_lat'): return 'ankle_pf_knee_flex'
    if has('soleus', 'tib_post', 'flex_dig', 'flex_hal', 'per_long', 'per_brev', 'peron'): return 'ankle_pf'
    if has('tib_ant', 'ext_dig', 'ext_hal', 'per_terr', 'per_tert'): return 'ankle_df'
    if has('ercspn', 'intobl', 'extobl'): return 'trunk'
    return 'other'

def build_synergy6():
    import numpy as np
    z = np.load(os.path.join(HERE, 'synergy_basis.npz'), allow_pickle=True)
    N, E = [], []
    def nd(t, l, x, y):
        N.append({'type': t, 'label': l, 'x': x, 'y': y}); return l
    def e(f, t, s, g, tag):
        E.append({'from': f, 'to': t, 'sign': s, 'gain': round(float(g), 3),
                  'tag': tag})
    nd('PORT-load', 'DRIVE / RG phase', 60, 300)
    nd('PORT-load', 'W fitted (fsa_backsolve rank-6)', 60, 380)
    for i, S in enumerate(('R', 'L')):
        X = 260 + i * 640
        names = [str(x) for x in z['muscle_names_' + S.lower()]]
        W = z['W_' + S.lower()]
        pools = {}
        for j, mn in enumerate(names):
            pools.setdefault(pool_of(mn), []).append((mn, W[j]))
        for k in range(6):
            syn = nd('IN-PF', 'Syn%d_%s' % (k + 1, S), X, 60 + k * 85)
            e('DRIVE / RG phase', syn, 'exc', 1.0, 'h(t) from H matrix')
            col = W[:, k]
            mx = col.max()
            shown = []
            for pname, members in pools.items():
                w = max(abs(m[1][k]) for m in members)
                if w >= 0.25 * mx and w > 0:
                    shown.append((pname, w, sum(m[1][k] for m in members)))
            for pname, w, _ in sorted(shown, key=lambda t: -t[1])[:5]:
                lbl = pname + '_' + S
                if not any(n['label'] == lbl for n in N):
                    existing = [n for n in N if n['label'].endswith('_' + S)
                                and n['type'] == 'MN']
                    y = 620 + (len(existing) % 8) * 62
                    nd('MN', lbl, X + 300, y)
                e(syn, lbl, 'exc', w, 'w=%.2f (pool max)' % w)
        for j in range(8):
            pass
    return {'nodes': N, 'edges': E,
            '_note': 'Rank-6 NMF synergies per leg (synergy_basis.npz, '
                     'VAF 0.945). Each Syn node drives muscle POOLS (max '
                     'member weight); full 43-muscle W matrix lives in the npz.'}


# ---------------------------------------------------------------- rules
def build_rules():
    """The 11 CONNECTOME.md rules as reference motifs."""
    N, E = [], []
    def nd(t, l, x, y):
        N.append({'type': t, 'label': l, 'x': x, 'y': y}); return l
    def e(f, t, s, g, tag):
        E.append({'from': f, 'to': t, 'sign': s, 'gain': g, 'tag': tag})
    R = [
        ('ia_homo', 'SN-Ia', 'Ia', 'IN-IaIN', None, 'MN', 'MN agonist', 'ia_to_mn'),
        ('ia_recip', 'SN-Ia', 'Ia', 'IN-IaIN', 'IaIN', 'MN', 'MN antagonist', 'ia_to_antagonist'),
        ('ii_exc', 'SN-II', 'II', 'IN-IIe', 'IIX', 'MN', 'MN agonist', 'ii_to_mn'),
        ('ii_inh', 'SN-II', 'II b', 'IN-IIi', 'IIIN', 'MN', 'MN antagonist', 'ia_to_antagonist'),
        ('ib_auto', 'SN-Ib', 'Ib', 'IN-IbIN', 'IBIN', 'MN', 'MN same', 'ib_to_mn_inh'),
        ('ib_rev', 'SN-Ib', 'Ib+', 'IN-Ib+', 'IBEXC', 'MN', 'MN ext (stance)', 'ib_group_exc'),
        ('heel', 'SN-heel', 'heel', 'IN-C', 'heel IN', 'HC-RG-E', 'RG-E (via InE/F)', 'heel_rge'),
        ('toe', 'SN-toe', 'toe', 'IN-C', 'toe IN', 'HC-RG-E', 'RG-E (via InE/F)', 'toe_rge'),
        ('ib_load', 'SN-Ib', 'Ib grp', 'IN-LBIN', 'LBIN', 'HC-RG-E', 'RG-E', 'ib_rge'),
        ('rc', 'MN', 'MN', 'RC', 'RC', 'MN', 'MN own+antag', 'renshaw'),
        ('cross', 'HC-RG-F', 'RG-F', 'IN-C', 'c1/V3', 'HC-RG-F', 'contra RG', 'c1_gain'),
    ]
    for i, (rid, st, sl, it, il, tt, tl, gk) in enumerate(R):
        x = 70 + (i % 4) * 300
        y = 90 + (i // 4) * 260
        s = nd(st, rid + ': ' + sl, x, y + 60)
        t = nd(tt, rid + ': ' + tl, x + 210, y + 60)
        e(s, t, 'exc' if rid not in ('ia_recip', 'ii_inh', 'ib_auto', 'rc') else 'inh',
          0.5, rid + ' (0 hops shown)')
        if il:
            m = nd(it, rid + ': ' + il, x + 105, y)
            e(s, m, 'exc', 1.0, rid)
            e(m, t, 'inh' if rid in ('ia_recip', 'ii_inh', 'ib_auto', 'rc')
              else 'exc', 0.5, rid + ' hops=1')
        # annotation node
        nd('PORT-load', rid + ': ' + gk, x + 105, y + 130)
    return {'nodes': N, 'edges': E,
            '_note': 'CONNECTOME.md rules: each motif = rule; tag carries the '
                     'gain key; hops=0 direct edge + hops=1 IN-layer edge both '
                     'shown where applicable. Ben edits via connectome_editor.'}


# ---------------------------------------------------------------- validate
def validate(name, spec):
    errs = []
    labels = set()
    for n in spec['nodes']:
        if n['type'] not in TYPES:
            errs.append('%s: bad type %s' % (name, n['type']))
        if n['label'] in labels:
            errs.append('%s: dup label %s' % (name, n['label']))
        labels.add(n['label'])
    for e in spec.get('edges', []):
        if e['from'] not in labels or e['to'] not in labels:
            errs.append('%s: dangling %s->%s' % (name, e['from'], e['to']))
    return errs


def main():
    out = {}
    errs = []
    out['biped2xcpg'] = build_biped2xcpg()
    out['biped2xcpg']['_note'] = (
        'source: Biped_2xCPG_wSubs\\connectivity_audit.txt (gains only, no '
        'polarity). Signs assigned by wiring-class convention (lamination '
        'INs/IaIN/RC inhibitory, afferent/PF/MN excitatory); edges not '
        'derivable that way are individually tagged "sign_assumed" — '
        'confirm against the aproj before using as ground truth. Muscle '
        'nodes + mn_to_muscle edges are IMPLIED (audit is neuron-only).')
    # W2L: prefer subagent mine, fall back to the coarse draft
    w2l = None
    for cand in (os.path.join(TMP, 'tpl_w2l.json'),
                 os.path.join(HERE, 'w2l_equivalent_draft.json')):
        if os.path.exists(cand):
            w2l = json.load(open(cand, encoding='utf-8'))
            if 'edges' not in w2l and 'synapses' in w2l:
                w2l['edges'] = w2l.pop('synapses')
            resolve_alias(w2l)
            out['w2laproj'] = w2l
            out['w2laproj']['_note'] = (
                'source: ' + os.path.basename(cand) +
                ' — 2023 ORIGINAL W2L: single L RG (endogenous, kickoff '
                'stimulus), NOT the contact-driven branch variant. Ia-chain '
                'typing by role: "Ia IN" = spindle front-end (SN-Ia), '
                '"Ia" = reciprocal-IN (IN-IaIN). All chemical gains 0.5; '
                'signs from SynapseType EquilibriumPotential.')
            break
    if w2l:
        out['bilateralrg'] = build_bilateralrg(w2l)
        out['bilateralrg']['_note'] = (
            'W2L base + documented BilateralRG chain (build_rg/comm/aff/'
            'contact .pl, tools/SESSION_NOTES_20260916.md): crossed L-RG->'
            'R-PF removed, R half-center added ipsilaterally + Stimulus_2 '
            'antiphase kickoff, Shinohara c1 (2.749) / V3 (0.1) commissurals, '
            'afferent-HC-excite 0.01, heel/toe contact C=20.')
    else:
        print('WARN: no w2l source yet -> bilateralrg skipped')
    li = os.path.join(TMP, 'tpl_li.json')
    if os.path.exists(li):
        out['li'] = json.load(open(li, encoding='utf-8'))
        out['li']['_note'] = (
            'source: Li Model\\walk tester rearranged.aproj mine '
            '(2026-09-23, spiking I&F generation). Contact-driven 2-layer '
            'CPG: per-leg stance CPG (6 nA tonic, tau_Ca 7.5 s burst '
            'termination) gated by CPG-inhibit; interleg = reciprocal '
            'stance-CPG inhibition + contralateral heel->CPG-inhibit + '
            'stance CPG->contra hip-swing PF. NO Renshaw/Ia/Ib/II chains '
            '(contact + hip-angle receptors only). Signs from SynapseType '
            'EquilibriumPotential; gains = per-synapse conductance (1-8 '
            'uS, not uniform).')
    else:
        print('WARN: tpl_li.json not ready yet -> li skipped (rerun later)')

    def jload(p):
        return json.load(open(os.path.join(HERE, p), encoding='utf-8'))
    v9 = jload('best_walk_params_v9.json')
    v10 = jload('best_walk_params_v10.json')
    s1 = jload('curriculum_stage1.json')
    s2 = jload('curriculum_stage2.json')
    s3 = jload('curriculum_stage3.json')
    s3k = jload(os.path.join('reports_20260923', 's3k_trial34_full_params.json'))
    for key, d, note in [
        ('walker_v9', v9, 'ground_walk_v10-era winner (study ground_walk_v9; pre-curriculum)'),
        ('walker_v10', v10, 'transient-reset + ankle-trim winner (ground_walk_v10_transient trial 56)'),
        ('walker_s1', s1, 'curriculum stage 1: air, deafferented (minimal knobs)'),
        ('walker_s2', s2, 'curriculum stage 2: air + afferents (heel/toe/central knobs emerged 0)'),
        ('walker_s3', s3, 'curriculum stage 3: ground (current pre-s3k family)'),
        ('walker_s3k', s3k, 'CURRENT production config (s3k study, trial 34; +pm phase machine)'),
    ]:
        out[key] = build_walker(key, d['params'],
                                pf_gain=d.get('pf_gain'), note=note)

    out['synergy6'] = build_synergy6()
    out['rules'] = build_rules()
    for nm in ('shevtsova', 'shinohara', 'rybak'):
        p = os.path.join(HERE, 'replication', nm + '_rules.json')
        out[nm] = convert_replication(p)
        out[nm]['_note'] = ('replication draft (Ben edits; ' + nm +
                            '_rules.json); weights = paper Table values')

    for nm, spec in out.items():
        errs += validate(nm, spec)
    if errs:
        print('VALIDATION ERRORS:')
        for x in errs[:40]:
            print(' ', x)
        sys.exit(1)
    dest = os.path.join(HERE, 'connectome_templates.json')
    json.dump(out, open(dest, 'w', encoding='utf-8'), indent=1)
    print('wrote', dest)
    for nm, spec in out.items():
        print('  %-14s %3d nodes %3d edges' %
              (nm, len(spec['nodes']), len(spec.get('edges', []))))


if __name__ == '__main__':
    main()
