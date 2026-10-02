"""Build the single-file SADb explorer app: app/sadb_app.html.

Reads export/sadb_export.json + export/sadb_cites.json + export/cluster_labels.json
and embeds everything into ONE self-contained HTML file (double-click to open;
core browsing works fully offline).

Tabs:
  Table   search/sort/filters + drill-through with a back-to-pivot bar
  Pivot   count matrix, labeled topic clusters, count/A-Z sortable rows+cols
  Map     year×citations + topic landscape; neuron styling; Research-Rabbit-style
          FOCUS MODE: click a paper -> BFS neighborhood to 1/2/3 degrees,
          DIRECT connections colored, 2nd-degree grayed; left panel = focus
          paper, right panel = connection list (click to hop)
  Search  ONLINE tab (new, 2026-09-28): live OpenAlex query + launch buttons for
          Google Scholar / PubMed / Web of Science (PSU) / doi.org — the app
          itself only fetches OpenAlex (CORS-open, no key); Scholar/WoS open in
          your browser so institutional access + their ToS apply

Re-run after any curation batch:
  myo python SADb_audit/export_corpus.py && SADb_audit/build_citation_graph.py && SADb_audit/app/build_app.py
"""
import datetime, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
records = json.load(open(os.path.join(SAD, "export", "sadb_export.json"), encoding="utf-8"))
cites = json.load(open(os.path.join(SAD, "export", "sadb_cites.json"), encoding="utf-8"))
cl_labels = {}
_lbl = os.path.join(SAD, "export", "cluster_labels.json")
if os.path.exists(_lbl):
    cl_labels = {int(k): v for k, v in json.load(open(_lbl, encoding="utf-8")).items()}

# --- source classification (same rule as before; export json lacks src) ---
import csv as _csv


def _rec_ids(path):
    ids = set()
    with open(path, encoding="utf-8-sig") as f:
        for row in _csv.reader(f):
            ids.update(c for c in row if c.startswith("rec") and len(c) == 17)
    return ids


_ORIG = _rec_ids(os.path.join(SAD, "airtable_papers_slim.csv"))
_DEMO = _rec_ids(os.path.join(SAD, "airtable_created50_ids.csv"))
_DIGEST = {"rechmhCcKA0ILQWKG", "recw5L1NKiDuwcdQY", "recp5dgFg2CMMDQAz",
           "rec8MVmUOx1oDk02j", "recB7n3thv2EEgswz", "rec1Q4YaKntrUybgY",
           "recCjTyA6Mz9tHw0G", "recZJU3Tv1NeOUDiG", "recEBU9QYruCYKihU",
           "recIjP2Mi5teZyaov", "recMZoiH13ytYevS7", "recp7yrg2a8E7DgyR",
           "recRmBChe7Lolwn0F", "recHoAYp5AUdgYEZA", "reck2RL4OlhuakV9q",
           "recafnPacssmJFhBt", "recenxdRrnk3ydqyg", "recoNlbHxaPTmt2Rs",
           "recpHeyTiV5C1veLV", "recl7iChHRGAfNn1c"}
_REST = set()
with open(os.path.join(SAD, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for _row in _csv.DictReader(f):
        if _row.get("DOI"):
            _REST.add(_row["DOI"].strip().lower())


def _classify(r):
    rid = r["id"]
    if rid in _ORIG:
        return "originals99"
    if rid in _DEMO:
        return "demo50"
    if rid in _DIGEST:
        return "digest20"
    if (r["doi"] or "").strip().lower() in _REST:
        return "rest383"
    return "task5import"


data = []
for r in records:
    data.append({
        "id": r["id"], "t": r["title"], "a": r["primary"], "s": "; ".join(r["secondary"]),
        "y": r["year"], "d": r["doi"], "an": r["animals"], "fb": r["feedback"],
        "af": r.get("afferents", []), "ast": r.get("animal_study", ""),
        "rs": r.get("robot_sim", ""), "pr": r.get("prune", ""),
        "md": r["models2"], "mr": r["models_ref"], "rv": r["reviews"],
        "pdf": 1 if r["has_pdf"] else 0, "n": 1 if r["has_notes"] else 0,
        "ar": 1 if r.get("archived") else 0,
        "nt": r["notes"], "src": _classify(r), "c": 0, "cl": -1,
        "lx": 0.0, "ly": 0.0,
    })

layout_path = os.path.join(SAD, "export", "sadb_layout.json")
if os.path.exists(layout_path):
    _lay = json.load(open(layout_path, encoding="utf-8"))
    for r in data:
        L = _lay.get(r["id"], {})
        r["c"] = L.get("c", 0)
        r["cl"] = L.get("cl", -1)
        r["lx"] = L.get("x", 0.0)
        r["ly"] = L.get("y", 0.0)
else:
    print("WARNING: sadb_layout.json missing — citation counts/clusters/landscape coords will be 0")

payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
cites_payload = json.dumps(cites, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
labels_payload = json.dumps(cl_labels, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SADb Explorer — Sensory Afferent Database</title>
<style>
 :root { --bg:#fafafa; --fg:#1a1a1a; --mut:#666; --line:#ddd; --acc:#0b62a4; }
 * { box-sizing:border-box; }
 body { font:14px/1.45 "Segoe UI",system-ui,sans-serif; margin:0; color:var(--fg); background:var(--bg); }
 header { padding:10px 16px; border-bottom:2px solid var(--acc); background:#fff; }
 h1 { font-size:17px; margin:0 0 2px; }
 .sub { color:var(--mut); font-size:12px; }
 .tabs { margin-top:8px; }
 .tabs button { border:1px solid var(--line); background:#fff; padding:6px 18px; cursor:pointer; font-size:14px; }
 .tabs button.on { background:var(--acc); color:#fff; border-color:var(--acc); }
 .bar { display:flex; flex-wrap:wrap; gap:8px; padding:8px 16px; background:#fff; border-bottom:1px solid var(--line); align-items:center; }
 .bar input[type=text], .bar select { padding:4px 8px; border:1px solid var(--line); }
 .bar input[type=number] { width:64px; padding:4px; border:1px solid var(--line); }
 .bar label { font-size:12px; color:var(--mut); white-space:nowrap; }
 main { padding:12px 16px; }
 table { border-collapse:collapse; width:100%; background:#fff; }
 th, td { border-bottom:1px solid var(--line); padding:5px 8px; text-align:left; vertical-align:top; }
 th { position:sticky; top:0; background:#eef3f8; cursor:pointer; user-select:none; white-space:nowrap; }
 td.num { text-align:right; white-space:nowrap; }
 tr.row:hover { background:#f0f6ff; cursor:pointer; }
 .pill { display:inline-block; background:#e7eef5; border-radius:8px; padding:0 7px; margin:1px 2px; font-size:11px; white-space:nowrap; }
 .pill.fb { background:#e8f2e4; } .pill.no { background:#fdeaea; } .pill.yes { background:#e4f0e8; }
 .pill.af { background:#f3e8f5; }
 #count { color:var(--mut); font-size:12px; margin:6px 0; }
 .pivot td, .pivot th { text-align:right; }
 .pivot td.rh, .pivot th.rh { text-align:left; }
 .pivot td.cell { cursor:pointer; }
 .pivot td.cell:hover { background:#ffe9b8; }
 .pivot td.tot { font-weight:600; background:#f4f4f4; }
 .phead { cursor:pointer; user-select:none; color:var(--acc); }
 #crumb { display:flex; align-items:center; gap:10px; padding:4px 0; font-size:13px; }
 #crumb button, .backbtn { border:1px solid var(--acc); background:#fff; color:var(--acc); padding:2px 10px; cursor:pointer; }
 #canvas { border:1px solid var(--line); background:#fff; cursor:grab; display:block; outline-offset:-3px; }
 #canvas:focus { outline:3px solid var(--acc); }
 #tip { position:fixed; display:none; background:#fff; border:1px solid var(--acc); padding:6px 9px;
        font-size:12px; max-width:340px; pointer-events:none; box-shadow:2px 2px 8px rgba(0,0,0,.15); z-index:50; }
 #detail { position:fixed; right:0; top:0; bottom:0; width:430px; background:#fff; border-left:2px solid var(--acc);
           padding:14px 16px; overflow-y:auto; display:none; z-index:60; box-shadow:-4px 0 16px rgba(0,0,0,.12); }
 #detail h2 { font-size:15px; margin:0 0 6px; } #detail .close { float:right; cursor:pointer; border:1px solid var(--line); padding:1px 8px; }
 #detail dt { font-weight:600; font-size:12px; color:var(--mut); margin-top:8px; }
 #detail dd { margin:2px 0 0; }
 .legend { display:flex; flex-wrap:wrap; gap:4px 12px; font-size:11px; color:var(--mut); margin-top:6px; }
 .legend i { display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:3px; }
 .warn { background:#fff4d5; border:1px solid #e8d49a; padding:4px 10px; font-size:12px; display:none; }
 .sr { position:absolute; width:1px; height:1px; overflow:hidden; clip:rect(0 0 0 0); }
 /* Research-Rabbit focus panels (map) */
 .fpanel { position:absolute; top:8px; bottom:8px; width:300px; background:rgba(255,255,255,.96);
           border:1px solid var(--acc); box-shadow:0 0 14px rgba(0,0,0,.15); overflow-y:auto;
           padding:10px 12px; font-size:12.5px; z-index:30; display:none; }
 #lpanel { left:8px; } #rpanel { right:8px; }
 .fpanel h3 { font-size:13px; margin:2px 0 6px; }
 .fpanel .ph { border-top:1px solid var(--line); margin-top:8px; padding-top:6px; font-weight:600; font-size:12px; color:var(--mut); }
 .conn { padding:3px 4px; border-radius:4px; cursor:pointer; }
 .conn:hover { background:#eef4fb; }
 .conn .ct { color:var(--mut); font-size:11px; }
 .conn.cur { background:#e2effb; }
 #view-map { position:relative; }
 .srch-r { border:1px solid var(--line); background:#fff; padding:8px 10px; margin:8px 0; border-radius:6px; }
 .srch-r .t { font-weight:600; }
 .srch-r .m { color:var(--mut); font-size:12px; margin:2px 0; }
 .srch-r a { color:var(--acc); margin-right:10px; font-size:12px; }
 .srch-r .incorpus { background:#e4f0e8; border-radius:8px; padding:0 7px; font-size:11px; }
 details summary { cursor:pointer; color:var(--acc); }
</style>
</head>
<body>
<header>
  <h1>SADb Explorer — Sensory Afferent Database <span class="sub" id="stats"></span></h1>
  <div class="sub">Snapshot of the Airtable "Sensory Feedback" corpus. Offline core; the Search tab adds live OpenAlex queries. Rebuild: <code>myo python SADb_audit/export_corpus.py &amp;&amp; SADb_audit/build_citation_graph.py &amp;&amp; SADb_audit/app/build_app.py</code></div>
  <div class="tabs">
    <button id="tab-table" class="on" onclick="show('table')">Table</button>
    <button id="tab-pivot" onclick="show('pivot')">Pivot</button>
    <button id="tab-map" onclick="show('map')">Bubble map</button>
    <button id="tab-search" onclick="show('search')">Search (online)</button>
  </div>
</header>

<div class="bar" id="bar-table">
  <input type="text" id="q" placeholder="Search title / author / year / DOI / notes…" size="38">
  <label>Animal <select id="f-an"><option value="">(any)</option></select></label>
  <label>Pathway <select id="f-fb"><option value="">(any)</option></select></label>
  <label>Afferent <select id="f-af"><option value="">(any)</option></select></label>
  <label>Source <select id="f-src"><option value="">(any)</option></select></label>
  <label>Year <input type="number" id="f-y0" placeholder="from"> – <input type="number" id="f-y1" placeholder="to"></label>
  <label><input type="checkbox" id="f-notes"> has notes</label>
  <label><input type="checkbox" id="f-pdf"> has PDF</label>
  <span class="warn" id="cellwarn"></span>
  <button class="backbtn" id="back-pivot" style="display:none" onclick="backToPivot()">◀ Back to pivot</button>
  <button onclick="clearCell()">Clear drill-through</button>
</div>

<div class="bar" id="bar-pivot" style="display:none">
  <label>Rows <select id="p-row"></select></label>
  <label>Columns <select id="p-col"></select></label>
  <label>Sort rows <select id="p-sr"><option value="count">by count</option><option value="alpha">A → Z</option><option value="alpha-r">Z → A</option></select></label>
  <label>Sort cols <select id="p-sc"><option value="count">by count</option><option value="alpha">A → Z</option><option value="alpha-r">Z → A</option></select></label>
  <label class="sub">Count = papers; multi-valued fields count once per value. Click a cell to drill through; the Table view then shows a "Back to pivot" button.</label>
</div>

<div class="bar" id="bar-map" style="display:none">
  <label>Layout <select id="m-layout">
    <option value="year">Year × citations</option>
    <option value="land">Topic landscape (citation network)</option>
  </select></label>
  <label>Style <select id="m-style">
    <option value="neurons">Neurons</option>
    <option value="plain">Plain bubbles</option>
  </select></label>
  <label>View <select id="m-view">
    <option value="all">All papers</option>
    <option value="refs">References</option>
    <option value="citedby">Cited by</option>
    <option value="review">Research review</option>
    <option value="animal">Animal studies</option>
    <option value="models">Models</option>
  </select></label>
  <label>Degrees <select id="m-deg">
    <option value="1">1 (direct only)</option>
    <option value="2" selected>2</option>
    <option value="3">3</option>
  </select></label>
  <label>Left click action <select id="m-swap">
    <option value="focus">Focus paper</option>
    <option value="details">Show details</option>
  </select></label>
  <label class="sub">Click = focus: network out to N degrees — <b>direct connections colored</b> (green triangles = papers citing the focus arrive excitatory; red circles = papers the focus cites), 2nd degree grayed. Left panel = focus paper; right panel = connections (click to hop). Esc = back.</label>
</div>

<div class="bar" id="bar-search" style="display:none">
  <input type="text" id="s-q" placeholder="Search OpenAlex live — title / author / topic…" size="44">
  <button onclick="doSearch()">Search OpenAlex</button>
  <label>per page <select id="s-n"><option>10</option><option selected>25</option><option>50</option></select></label>
  <span class="sub" id="s-note"></span>
</div>

<main id="view-table">
  <div id="count"></div>
  <table><thead><tr>
    <th data-k="a">First author</th><th data-k="t">Title</th><th data-k="y" class="num">Year</th>
    <th data-k="c" class="num">Cites</th><th>Tags</th>
  </tr></thead><tbody id="tbody"></tbody></table>
</main>

<main id="view-pivot" style="display:none"></main>

<main id="view-map" style="display:none">
  <div id="crumb"></div>
  <canvas id="canvas" tabindex="0" role="application"
          aria-label="Bubble map. Use arrow keys to move between bubbles, Enter for the primary action, Shift+Enter for details, Escape to go back."></canvas>
  <div id="lpanel" class="fpanel" aria-label="Focus paper details"></div>
  <div id="rpanel" class="fpanel" aria-label="Connection list"></div>
  <div class="sr" id="map-live" aria-live="polite"></div>
  <div class="legend" id="maplegend"></div>
</main>

<main id="view-search" style="display:none">
  <div class="sub" style="margin:4px 0 10px">
    Live literature search against <b>OpenAlex</b> (open, no key). For sources that need your
    institutional session or have terms of service, each result has launch buttons that open the
    query in your own browser: <b>Google Scholar</b>, <b>PubMed</b>, and <b>Web of Science</b>
    (opens the PSU-libraries search page with the query copied to your clipboard).
    Corpus DOIs already in this snapshot are badged. Works offline for everything except the
    OpenAlex fetch itself.
  </div>
  <div id="s-results"><div class="sub">Type a query and press Search.</div></div>
</main>

<div id="tip"></div>
<div id="detail"><span class="close" onclick="hideDetail()">✕</span><div id="detail-body"></div></div>

<script type="application/json" id="sadb-data">__DATA__</script>
<script type="application/json" id="sadb-cites">__CITES__</script>
<script type="application/json" id="cluster-labels">__CLUSTER__</script>
<script>
const DATA = JSON.parse(document.getElementById('sadb-data').textContent);
const CITES = JSON.parse(document.getElementById('sadb-cites').textContent);
const CL_LABELS = JSON.parse(document.getElementById('cluster-labels').textContent);
/* Ben's accessible palette from Code\Matlab\Colors.m */
const PAL = ['#FFD700','#FFB14E','#FA8775','#EA5F94','#CD34B5','#9D02D7','#0000FF'];
const GREY = '#c9c9c9';
const SRC_LABEL = {originals99:'original 99', demo50:'demo 50', rest383:'rest import 383',
                   digest20:'RR digest 20', task5import:'task5 import'};
const byId = {}; DATA.forEach(r => byId[r.id] = r);
const citedBy = {};
for (const [p, refs] of Object.entries(CITES)) refs.forEach(t => { (citedBy[t] = citedBy[t] || []).push(p); });
const neighborsOf = id => [...new Set([...(CITES[id]||[]), ...(citedBy[id]||[])])];
const ANIMALS = [...new Set(DATA.flatMap(r => r.an))].sort();
const AFFECTS = [...new Set(DATA.flatMap(r => r.af))].sort();
function clName(cl){ return CL_LABELS[String(cl)] || CL_LABELS[cl] || ('Cluster '+(cl+1)); }
function clShort(cl){ const s = clName(cl); return s.length>34 ? s.slice(0,32)+'…' : s; }

const state = { q:'', an:'', fb:'', af:'', src:'', y0:null, y1:null, notes:false, pdf:false,
                sort:{k:'c', dir:-1}, cell:null };
const $ = id => document.getElementById(id);
function esc(s){ return (s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;'); }

// ================= TABLE =================
function hay(r){ return (r.t+' '+r.a+' '+r.s+' '+r.y+' '+r.d+' '+r.nt).toLowerCase(); }
function matches(r){
  const st = state;
  if (st.q && !hay(r).includes(st.q.toLowerCase())) return false;
  if (st.an && !r.an.includes(st.an)) return false;
  if (st.fb && !r.fb.includes(st.fb)) return false;
  if (st.af && !r.af.includes(st.af)) return false;
  if (st.src && r.src !== st.src) return false;
  if (st.y0 && !(r.y >= st.y0)) return false;
  if (st.y1 && !(r.y <= st.y1)) return false;
  if (st.notes && !r.n) return false;
  if (st.pdf && !r.pdf) return false;
  if (st.cell){
    const [rk, rv, ck, cv] = st.cell;
    if (!dimHas(r, rk, rv) || !dimHas(r, ck, cv)) return false;
  }
  return true;
}
function filtered(){
  const out = DATA.filter(matches);
  const {k, dir} = state.sort;
  out.sort((a,b) => {
    let va = a[k], vb = b[k];
    if (k === 'y'){ va = va||0; vb = vb||0; }
    return (va>vb?1:va<vb?-1:0)*dir;
  });
  return out;
}

const DIMS = {
  cluster:  { label:'Topic cluster', vals:r=>[r.cl>=0? clName(r.cl) : 'no citation links'] },
  decade:  { label:'Decade', vals:r=>[r.y? Math.floor(r.y/10)*10+'s' : 'no year'] },
  author:  { label:'Primary author', vals:r=>[r.a||'(none)'] },
  animal:  { label:'Animal', multi:true, vals:r=> r.an.length? r.an : ['(untagged)'] },
  pathway: { label:'Feedback pathway', multi:true, vals:r=> r.fb.length? r.fb : ['(no pathway link)'] },
  afferent:{ label:'Afferent types', multi:true, vals:r=> r.af.length? r.af : ['(untagged)'] },
  notes:   { label:'Notes', vals:r=>[r.n? 'has notes':'no notes'] },
  pdf:     { label:'PDF', vals:r=>[r.pdf? 'has PDF':'no PDF'] },
  source:  { label:'Import source (legacy)', vals:r=>[SRC_LABEL[r.src]||r.src||'?'] },
};
function dimHas(r, key, val){ return DIMS[key].vals(r).includes(val); }

function renderTable(){
  const rows = filtered();
  $('count').textContent = rows.length + ' of ' + DATA.length + ' papers' +
    (state.cell ? '  [drill-through active]' : '');
  $('back-pivot').style.display = state.cell ? 'inline-block' : 'none';
  const tb = $('tbody'); tb.innerHTML = '';
  const frag = document.createDocumentFragment();
  rows.slice(0, 400).forEach(r => {
    const tr = document.createElement('tr'); tr.className = 'row';
    const tags = [];
    if (r.n) tags.push('<span class="pill yes">notes</span>');
    else tags.push('<span class="pill no">no notes</span>');
    if (r.ar) tags.push('<span class="pill" style="background:#eee3f3">app-only (removed from Airtable)</span>');
    if (r.pdf) tags.push('<span class="pill yes">PDF</span>');
    if (r.pr === 'prune-candidate') tags.push('<span class="pill no">prune?</span>');
    r.an.slice(0,3).forEach(a => tags.push('<span class="pill">'+a+'</span>'));
    r.af.slice(0,3).forEach(a => tags.push('<span class="pill af">'+esc(a)+'</span>'));
    r.fb.slice(0,2).forEach(f => tags.push('<span class="pill fb">'+f+'</span>'));
    tr.innerHTML = '<td>'+(r.a||'')+'</td><td>'+esc(r.t)+
      '<div class="sub">'+esc(r.s)+'</div></td><td class="num">'+(r.y||'')+'</td>'+
      '<td class="num">'+r.c+'</td><td>'+tags.join(' ')+'</td>';
    tr.onclick = () => showDetail(r);
    frag.appendChild(tr);
  });
  tb.appendChild(frag);
  if (rows.length > 400) $('count').textContent += ' — showing first 400 (use filters)';
}
document.querySelectorAll('th[data-k]').forEach(th => th.onclick = () => {
  const k = th.dataset.k;
  state.sort = {k, dir: state.sort.k===k ? -state.sort.dir : 1};
  renderTable();
});
['q','f-an','f-fb','f-af','f-src','f-y0','f-y1','f-notes','f-pdf'].forEach(id => {
  $(id).addEventListener(id==='q' ? 'input' : 'change', () => {
    state.q = $('q').value; state.an = $('f-an').value; state.fb = $('f-fb').value;
    state.af = $('f-af').value; state.src = $('f-src').value;
    state.y0 = $('f-y0').value ? +$('f-y0').value : null;
    state.y1 = $('f-y1').value ? +$('f-y1').value : null;
    state.notes = $('f-notes').checked; state.pdf = $('f-pdf').checked;
    renderTable(); if (curTab==='map') drawMap();
  });
});
function clearCell(){ state.cell = null; $('cellwarn').style.display='none';
  $('back-pivot').style.display='none'; renderTable(); }
function backToPivot(){ clearCell(); show('pivot'); }

// ================= PIVOT =================
function sortEntries(entries, mode){
  if (mode === 'alpha') return entries.sort((a,b)=> a[0].localeCompare(b[0]));
  if (mode === 'alpha-r') return entries.sort((a,b)=> b[0].localeCompare(a[0]));
  return entries.sort((a,b)=> b[1]-a[1]);
}
function renderPivot(){
  const rk = $('p-row').value, ck = $('p-col').value;
  const rVals = new Map(), cVals = new Map(), matrix = new Map();
  DATA.forEach(r => {
    const rvs = DIMS[rk].vals(r), cvs = DIMS[ck].vals(r);
    rvs.forEach(rv => rVals.set(rv,(rVals.get(rv)||0)+1));
    cvs.forEach(cv => cVals.set(cv,(cVals.get(cv)||0)+1));
    rvs.forEach(rv => cvs.forEach(cv => {
      const key = rv+'||'+cv; matrix.set(key,(matrix.get(key)||0)+1);
    }));
  });
  const rl = sortEntries([...rVals.entries()], $('p-sr').value).slice(0,40);
  const cl = sortEntries([...cVals.entries()], $('p-sc').value).slice(0,24);
  let h = '<table class="pivot"><tr><th class="rh">'+DIMS[rk].label+' ↓ / '+DIMS[ck].label+' →</th>';
  cl.forEach(([cv]) => h += '<th>'+esc(cv)+'</th>');
  h += '<th>total</th></tr>';
  rl.forEach(([rv]) => {
    h += '<tr><td class="rh">'+esc(rv)+'</td>';
    let tot = 0;
    cl.forEach(([cv]) => {
      const v = matrix.get(rv+'||'+cv)||0; tot += v;
      h += '<td class="cell" data-r="'+esc(rv)+'" data-c="'+esc(cv)+'"'+(v?'':' style="color:#ccc"')+'>'+(v||'·')+'</td>';
    });
    h += '<td class="tot">'+tot+'</td></tr>';
  });
  h += '</table>';
  if (rVals.size>40 || cVals.size>24) h += '<p class="sub">Showing top 40 row / 24 column values ('+rVals.size+' × '+cVals.size+' total).</p>';
  const el = $('view-pivot'); el.innerHTML = h;
  el.querySelectorAll('td.cell').forEach(td => td.onclick = () => {
    state.cell = [rk, td.dataset.r, ck, td.dataset.c];
    $('cellwarn').textContent = 'Drill: '+td.dataset.r+' × '+td.dataset.c;
    $('cellwarn').style.display = 'inline';
    show('table');
  });
}

// ================= MAP =================
const cv = $('canvas'), ctx = cv.getContext('2d');
const MAXC = Math.max(...DATA.map(r=>r.c));
const PAD = {l:62, r:24, t:18, b:40};
const X0 = 1905, X1 = 2030, Y0 = 0, Y1 = 5.3;   // year, log10(citations+1)
let view = {k:1, tx:0, ty:0};
let drag = null, dragMoved = 0, hover = null, kfocus = -1;
const ms = { mode:'all', focus:null, spot:null };

function hashId(s){ let h = 2166136261; for (let i=0;i<s.length;i++){ h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); } return h>>>0; }
const POS = {}, RAD = {}, LPOS = {};
DATA.forEach(r => {
  const jh = hashId(r.id);
  POS[r.id] = [
    (r.y || 1970) + ((jh % 1000)/1000 - 0.5) * 1.8,
    Math.log10((r.c||0) + 1) + (((jh>>>10) % 1000)/1000 - 0.5) * 0.16,
  ];
  LPOS[r.id] = [r.lx, r.ly];
  RAD[r.id] = Math.min(24, 1.7 + 5.5*Math.pow((r.c||0)/MAXC, 0.25)*Math.pow(MAXC, 0.12));
});
function bx(x){ return PAD.l + (x - X0)/(X1 - X0) * (cv.width - PAD.l - PAD.r); }
function by(y){ return cv.height - PAD.b - (y - Y0)/(Y1 - Y0) * (cv.height - PAD.t - PAD.b); }
function radOf(r){ return RAD[r.id]; }

// focus-mode BFS: degree map over the undirected in-corpus citation graph
let degMap = null;   // {id: degree} when a focus neighborhood is active
function computeDegMap(fid, maxDeg){
  const dm = {[fid]: 0};
  let frontier = [fid];
  for (let d = 1; d <= maxDeg; d++){
    const next = [];
    frontier.forEach(id => neighborsOf(id).forEach(nb => {
      if (!(nb in dm)){ dm[nb] = d; next.push(nb); }
    }));
    frontier = next;
    if (Object.keys(dm).length > 800) break;   // hub guard
  }
  return dm;
}

// PX: canvas-pixel position per paper, per layout mode
let curLayout = 'year';
const PX = {};
function recomputePX(){
  if (curLayout === 'land'){
    let x0=1e9,x1=-1e9,y0=1e9,y1=-1e9;
    DATA.forEach(r => { const [a,b] = LPOS[r.id];
      if(a<x0)x0=a; if(a>x1)x1=a; if(b<y0)y0=b; if(b>y1)y1=b; });
    const sx = (x1-x0)||1, sy = (y1-y0)||1, m = 30;
    DATA.forEach(r => {
      PX[r.id] = [ m + (LPOS[r.id][0]-x0)/sx*(cv.width-2*m),
                   cv.height - m - (LPOS[r.id][1]-y0)/sy*(cv.height-2*m) ];
    });
  } else {
    DATA.forEach(r => { PX[r.id] = [bx(POS[r.id][0]), by(POS[r.id][1])]; });
  }
}

function nodesNow(){
  if (ms.spot && degMap){
    return DATA.filter(r => r.id in degMap);
  }
  if ((ms.mode==='refs' || ms.mode==='citedby') && ms.focus){
    const nb = ms.mode==='refs' ? (CITES[ms.focus]||[]) : (citedBy[ms.focus]||[]);
    const set = new Set([ms.focus, ...nb]);
    return DATA.filter(r => set.has(r.id));
  }
  return DATA;
}
function colorOf(r){
  if (ms.spot && degMap){
    const d = degMap[r.id];
    if (d >= 2) return GREY;                       // 2nd degree grayed (Ben)
    // direct connections keep their semantic color
  }
  const m = ms.mode;
  if (m==='review') return r.rv.length ? PAL[5] : GREY;
  if (m==='models') return r.md.length ? PAL[4] : (r.mr.length ? PAL[1] : GREY);
  if (m==='animal'){
    if (!r.an.length) return GREY;
    return PAL[ANIMALS.indexOf(r.an[0]) % PAL.length];
  }
  return r.cl >= 0 ? PAL[r.cl % PAL.length] : GREY;
}
function sizeCanvas(){
  cv.width = Math.min(1600, window.innerWidth-60); cv.height = Math.min(880, window.innerHeight-200);
}
function drawMap(){
  sizeCanvas();
  ctx.setTransform(1,0,0,1,0,0);
  ctx.fillStyle = '#fff'; ctx.fillRect(0,0,cv.width,cv.height);
  if (curLayout === 'year'){
    ctx.strokeStyle = '#e3e3e3'; ctx.fillStyle = '#777';
    ctx.font = '11px system-ui'; ctx.lineWidth = 1;
    for (let yr = 1920; yr <= 2020; yr += 20){
      const X = bx(yr);
      ctx.beginPath(); ctx.moveTo(X, PAD.t); ctx.lineTo(X, cv.height-PAD.b); ctx.stroke();
      ctx.textAlign = 'center'; ctx.fillText(String(yr), X, cv.height-PAD.b+16);
    }
    const ylab = ['1','10','100','1k','10k','100k'];
    for (let i = 0; i <= 5; i++){
      const Y = by(i);
      ctx.beginPath(); ctx.moveTo(PAD.l, Y); ctx.lineTo(cv.width-PAD.r, Y); ctx.stroke();
      ctx.textAlign = 'right'; ctx.fillText(ylab[i], PAD.l-8, Y+4);
    }
    ctx.textAlign = 'center'; ctx.fillStyle = '#444';
    ctx.fillText('Publication year', cv.width/2, cv.height-6);
    ctx.save(); ctx.translate(14, cv.height/2); ctx.rotate(-Math.PI/2);
    ctx.fillText('Citations (log scale)', 0, 0); ctx.restore();
  } else {
    ctx.fillStyle = '#888'; ctx.font = '12px system-ui'; ctx.textAlign = 'center';
    ctx.fillText('Topic landscape — proximity ≈ citation similarity ('+
      Object.keys(CL_LABELS).length+' named clusters)', cv.width/2, 14);
  }

  const shown = nodesNow();
  const searchSet = state.q ? new Set(filtered().map(r=>r.id)) : null;
  const ordered = [...shown].sort((a,b) => RAD[b.id]-RAD[a.id]);
  const dimmed = r => (searchSet && !searchSet.has(r.id));
  const neuronStyle = $('m-style').value === 'neurons';
  ctx.setTransform(view.k,0,0,view.k,view.tx,view.ty);
  if (ms.spot && degMap) drawNeighborhoodEdges();
  if (neuronStyle) ordered.forEach(r => { if (!dimmed(r)) drawNeuronParts(r); });
  ordered.forEach(r => {
    const [sx, sy] = PX[r.id];
    ctx.beginPath(); ctx.arc(sx, sy, RAD[r.id], 0, 6.2832);
    let alpha = 0.85;
    if (dimmed(r)) alpha = 0.10;
    if (ms.spot && degMap && degMap[r.id] >= 2) alpha = 0.38;   // grayed 2nd degree
    ctx.globalAlpha = alpha;
    ctx.fillStyle = colorOf(r); ctx.fill();
    if (r.id === (ms.spot||'')){ ctx.globalAlpha = 1; ctx.lineWidth = 3.5/view.k; ctx.strokeStyle = '#111'; ctx.stroke(); }
    if (r.id === (ms.focus||'')){ ctx.globalAlpha = 1; ctx.lineWidth = 2.5/view.k; ctx.strokeStyle = '#111'; ctx.stroke(); }
    if (r === hover || (kfocus >= 0 && r === kFocused())) {
      ctx.globalAlpha = 1; ctx.lineWidth = 2/view.k; ctx.strokeStyle = '#000'; ctx.setLineDash([4/view.k,3/view.k]); ctx.stroke(); ctx.setLineDash([]);
    }
  });
  if (neuronStyle) drawSynapses();
  ctx.globalAlpha = 1;
  ctx.setTransform(1,0,0,1,0,0);
  renderLegend(); updateCrumb(); renderPanels();
}

// faint undirected edges across the whole shown neighborhood (focus mode)
function drawNeighborhoodEdges(){
  ctx.lineWidth = 0.8/view.k; ctx.globalAlpha = 0.22; ctx.strokeStyle = '#9aa4ad';
  let drawn = 0;
  for (const [p, refs] of Object.entries(CITES)){
    if (!(p in degMap)) continue;
    for (const t of refs){
      if ((t in degMap) && t !== p){
        const A = PX[p], B = PX[t];
        if (A && B){ ctx.beginPath(); ctx.moveTo(A[0],A[1]); ctx.lineTo(B[0],B[1]); ctx.stroke();
          if (++drawn > 3000) { ctx.globalAlpha = 1; return; } }
      }
    }
  }
  ctx.globalAlpha = 1;
}

// --- neuron morphology (dendrites + axon stub), deterministic per paper ---
function drawNeuronParts(r){
  const [sx, sy] = PX[r.id], rad = RAD[r.id];
  const h = hashId(r.id), col = colorOf(r);
  ctx.strokeStyle = col; ctx.globalAlpha = 0.5; ctx.lineWidth = 1.0/view.k;
  const nd = 4 + (h % 3);
  for (let i = 0; i < nd; i++){
    const a = (h % 628)/100 + i * 6.2832/nd + ((h>>>(4+i))%100)/300;
    const L1 = rad*(1.4 + ((h>>>(8+i))%100)/100);
    const x1 = sx + Math.cos(a)*(rad+L1), y1 = sy + Math.sin(a)*(rad+L1);
    const a2 = a + (((h>>>(12+i))%100)/100 - 0.5);
    ctx.beginPath();
    ctx.moveTo(sx + Math.cos(a)*rad*0.9, sy + Math.sin(a)*rad*0.9);
    ctx.lineTo(x1, y1);
    ctx.lineTo(x1 + Math.cos(a2)*L1*0.4, y1 + Math.sin(a2)*L1*0.4);
    ctx.stroke();
  }
  const aa = ((h>>>20)%628)/100;
  const ax = sx + Math.cos(aa)*rad*2.6, ay = sy + Math.sin(aa)*rad*2.6;
  const mx = (sx+ax)/2 + Math.cos(aa+1.2)*rad*0.4, my = (sy+ay)/2 + Math.sin(aa+1.2)*rad*0.4;
  ctx.beginPath(); ctx.moveTo(sx + Math.cos(aa)*rad, sy + Math.sin(aa)*rad);
  ctx.quadraticCurveTo(mx, my, ax, ay); ctx.stroke();
  ctx.beginPath(); ctx.arc(ax, ay, Math.max(1.2/view.k, rad*0.18), 0, 6.2832); ctx.stroke();
  ctx.globalAlpha = 1;
}

// --- synapses of the focused/spotlighted paper (Ben's convention: open triangle
// = excitatory, filled circle = inhibitory) ---
function pathBetween(A, rA, B, rB){
  const d = Math.hypot(B[0]-A[0], B[1]-A[1]) || 1;
  const ux = (B[0]-A[0])/d, uy = (B[1]-A[1])/d;
  const px = -uy, py = ux;
  ctx.beginPath();
  ctx.moveTo(A[0]+ux*(rA+2), A[1]+uy*(rA+2));
  ctx.quadraticCurveTo((A[0]+B[0])/2 + px*d*0.10, (A[1]+B[1])/2 + py*d*0.10,
                       B[0]-ux*(rB+6), B[1]-uy*(rB+6));
  ctx.stroke();
}
function drawSynapses(){
  const fid = ms.spot || ms.focus; if (!fid || !PX[fid]) return;
  const F = PX[fid];
  const inc = (citedBy[fid]||[]).filter(id => PX[id] && (!degMap || degMap[id] === 1)).slice(0, 400);
  const out = (CITES[fid]||[]).filter(id => PX[id] && (!degMap || degMap[id] === 1)).slice(0, 400);
  ctx.lineWidth = 1.6/view.k; ctx.globalAlpha = 0.8;
  ctx.strokeStyle = '#2e7d32';
  inc.forEach(id => pathBetween(PX[id], RAD[id], F, RAD[fid]));
  ctx.strokeStyle = '#b03040';
  out.forEach(id => pathBetween(F, RAD[fid], PX[id], RAD[id]));
  const t = Math.max(3.5, 7/view.k);
  inc.forEach(id => triangleAt(F, PX[id], t));
  out.forEach(id => {
    const P = PX[id];
    const d = Math.hypot(F[0]-P[0], F[1]-P[1]) || 1;
    ctx.beginPath();
    ctx.arc(P[0] + (F[0]-P[0])/d*(RAD[id]+t*0.9), P[1] + (F[1]-P[1])/d*(RAD[id]+t*0.9),
            t*0.7, 0, 6.2832);
    ctx.fillStyle = '#b03040'; ctx.fill();
  });
  ctx.globalAlpha = 1;
}
function triangleAt(F, P, t){
  const d = Math.hypot(F[0]-P[0], F[1]-P[1]) || 1;
  const ux = (F[0]-P[0])/d, uy = (F[1]-P[1])/d;
  const bx = F[0] - ux*(t*1.6), by = F[1] - uy*(t*1.6);
  ctx.beginPath();
  ctx.moveTo(F[0] - ux*1.2, F[1] - uy*1.2);
  ctx.lineTo(bx - uy*t, by + ux*t);
  ctx.lineTo(bx + uy*t, by - ux*t);
  ctx.closePath();
  ctx.fillStyle = '#fff'; ctx.fill();
  ctx.strokeStyle = '#2e7d32'; ctx.lineWidth = 1.1/view.k; ctx.stroke();
}
function kFocused(){ return kfocus >= 0 ? nodesSorted()[kfocus] : null; }
function nodesSorted(){
  const ns = nodesNow();
  if (!nodesSorted._c || nodesSorted._c !== ns){
    nodesSorted._list = [...ns].sort((a,b)=> (a.y||0)-(b.y||0) || a.c-b.c);
    nodesSorted._c = ns;
  }
  return nodesSorted._list;
}

// --- Research-Rabbit panels ---
function connRow(id, mark){
  const r = byId[id]; if (!r) return '';
  return '<div class="conn'+(mark?' cur':'')+'" onclick="hopFocus(\''+id+'\')">'+
    '<b>'+esc(r.a||'?')+'</b> '+(r.y||'')+
    (mark? ' <span class="ct">['+mark+']</span>':'')+
    '<div class="ct">'+esc(r.t.slice(0,72))+(r.t.length>72?'…':'')+'</div></div>';
}
function renderPanels(){
  const L = $('lpanel'), R = $('rpanel');
  if (!(ms.spot && degMap)){ L.style.display='none'; R.style.display='none'; return; }
  const r = byId[ms.spot]; if (!r) return;
  L.style.display='block'; R.style.display='block';
  const nRef = (CITES[r.id]||[]).length, nBy = (citedBy[r.id]||[]).length;
  L.innerHTML = '<h3>Focus paper</h3>'+
    '<b>'+esc(r.t)+'</b><div class="ct">'+esc(r.a)+(r.s? ' — '+esc(r.s):'')+' · '+(r.y||'n.d.')+
    ' · '+r.c+' citations</div>'+
    (r.d? '<div style="margin-top:6px"><a href="https://doi.org/'+esc(r.d)+'" target="_blank" rel="noopener">'+esc(r.d)+'</a></div>':'')+
    '<div class="ph">In-corpus graph</div>'+
    '<div>'+nRef+' references · '+nBy+' citing papers · '+
      Object.keys(degMap).length+' shown ('+$('m-deg').value+' degrees)</div>'+
    '<div class="ph">Animals</div><div>'+(r.an.length? r.an.map(a=>'<span class="pill">'+esc(a)+'</span>').join(''):'—')+'</div>'+
    '<div class="ph">Afferents</div><div>'+(r.af.length? r.af.map(a=>'<span class="pill af">'+esc(a)+'</span>').join(''):'—')+'</div>'+
    '<div class="ph">Pathways</div><div>'+(r.fb.length? r.fb.map(f=>'<span class="pill fb">'+esc(f)+'</span>').join(''):'—')+'</div>'+
    '<div class="ph">Actions</div>'+
    '<button class="backbtn" onclick="viewRefs(\''+r.id+'\')">References map</button> '+
    '<button class="backbtn" onclick="viewCitedBy(\''+r.id+'\')">Cited-by map</button> '+
    '<button class="backbtn" onclick="showDetail(byId[\''+r.id+'\'])">Details</button> '+
    '<button class="backbtn" onclick="clearSpot()">✕ Un-focus</button>'+
    (r.n? '<div class="ph">Curation note</div><div>'+esc(r.nt.slice(0,420))+'</div>':'')+
    (r.rs? '<div class="ph">Robot/sim translation</div><div>'+esc(r.rs.slice(0,420))+'</div>':'');
  const cap = 60;
  const refs = (CITES[r.id]||[]), bys = (citedBy[r.id]||[]);
  const second = Object.keys(degMap).filter(id => degMap[id] >= 2);
  R.innerHTML = '<h3>Connections</h3>'+
    '<div class="ph">▼ Cited by '+bys.length+' (excitatory input)</div>'+
    (bys.length? bys.slice(0,cap).map(id => connRow(id)).join('') +
      (bys.length>cap? '<div class="ct">+'+(bys.length-cap)+' more…</div>':'') : '<div class="ct">none in corpus</div>')+
    '<div class="ph">▼ References '+refs.length+' (focus cites)</div>'+
    (refs.length? refs.slice(0,cap).map(id => connRow(id)).join('') +
      (refs.length>cap? '<div class="ct">+'+(refs.length-cap)+' more…</div>':'') : '<div class="ct">none in corpus</div>')+
    '<div class="ph">2nd degree ('+second.length+')</div><div class="ct">'+
      second.slice(0,40).map(id => (byId[id]? (byId[id].a||'?')+' '+(byId[id].y||''):'')).join(' · ')+
      (second.length>40? ' …':'')+'</div>';
}
function hopFocus(id){
  ms.spot = id; kfocus = -1; nodesSorted._c = null;
  degMap = computeDegMap(id, +$('m-deg').value);
  resetView();
}
function clearSpot(){ ms.spot = null; degMap = null; kfocus = -1; nodesSorted._c = null; resetView(); }

function renderLegend(){
  const m = ms.mode; let items = [];
  if (m==='review') items = [['has Review Papers link', PAL[5]], ['not tagged as review', GREY]];
  else if (m==='models') items = [['is a model (Is-the-model-paper link)', PAL[4]], ['paper references models', PAL[1]], ['neither', GREY]];
  else if (m==='animal'){
    const used = {}; DATA.forEach(r => { if (r.an.length) used[r.an[0]] = (used[r.an[0]]||0)+1; });
    items = Object.entries(used).sort((a,b)=>b[1]-a[1]).slice(0,14)
      .map(([a,n]) => [a+' ('+n+')', PAL[ANIMALS.indexOf(a)%PAL.length]]);
  } else {
    const used = {}; DATA.forEach(r => { if (r.cl>=0) used[r.cl] = (used[r.cl]||0)+1; });
    items = Object.keys(used).sort((a,b)=>used[b]-used[a]).slice(0,14)
      .map(c => [clShort(+c)+' ('+used[c]+')', PAL[c%PAL.length]]);
    if (ms.spot && degMap) items.push(['2nd degree (grayed)', GREY]);
  }
  $('maplegend').innerHTML = items.map(([t,c])=>'<span><i style="background:'+c+'"></i>'+esc(t)+'</span>').join('');
}
function updateCrumb(){
  const c = $('crumb');
  if ((ms.mode==='refs'||ms.mode==='citedby') && ms.focus){
    const r = byId[ms.focus];
    const n = ms.mode==='refs' ? (CITES[ms.focus]||[]).length : (citedBy[ms.focus]||[]).length;
    c.innerHTML = '<button id="map-back">← All papers</button> '+
      (ms.mode==='refs'?'References of':'Cited by')+' <b>'+esc((r.a||'')+' '+(r.y||''))+'</b>: '+
      esc(r.t.slice(0,80))+' — '+n+' in-corpus papers';
  } else if (ms.spot){
    const r = byId[ms.spot];
    c.innerHTML = '<button id="map-back">← Show all</button> Focused: <b>'+esc((r.a||'')+' '+(r.y||''))+'</b> — '+
      esc(r.t.slice(0,80))+' · '+Object.keys(degMap||{}).length+' papers within '+$('m-deg').value+' degree(s)';
  } else {
    const lbl = {all:'All papers', review:'Research review', animal:'Animal studies', models:'Models'}[ms.mode];
    c.innerHTML = '<span class="sub">'+lbl+' — click a bubble to focus its citation neighborhood; switch View to References or Cited by for one-direction maps.</span>';
  }
  const back = $('map-back');
  if (back) back.onclick = () => {
    if (ms.mode==='refs'||ms.mode==='citedby'){ ms.focus = null; } else { ms.spot = null; degMap = null; }
    kfocus = -1; nodesSorted._c = null; resetView();
  };
}
function resetView(){
  sizeCanvas();
  recomputePX();
  view = {k:1, tx:0, ty:0};
  const ns = nodesNow();
  const xs = ns.map(r=>PX[r.id][0]), ys = ns.map(r=>PX[r.id][1]);
  const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  view.k = Math.max(0.4, Math.min(cv.width/(x1-x0+80), cv.height/(y1-y0+80), 4));
  view.tx = cv.width/2 - (x0+x1)/2*view.k; view.ty = cv.height/2 - (y0+y1)/2*view.k;
  drawMap();
}
function primaryAct(r){
  if (ms.mode==='refs'){ ms.focus = r.id; ms.spot = null; degMap = null; kfocus = -1; nodesSorted._c = null; resetView(); return; }
  if (ms.mode==='citedby'){ ms.focus = r.id; ms.spot = null; degMap = null; kfocus = -1; nodesSorted._c = null; resetView(); return; }
  hopFocus(r.id);
}
function secondaryAct(r){ showDetail(r); }
function act(r, primary){
  const doFocus = $('m-swap').value==='focus' ? primary : !primary;
  if (doFocus) primaryAct(r); else secondaryAct(r);
}

function hitTest(mx, my){
  const ordered = [...nodesNow()].sort((a,b)=>RAD[a.id]-RAD[b.id]);
  for (const r of ordered){
    const sx = PX[r.id][0]*view.k + view.tx, sy = PX[r.id][1]*view.k + view.ty;
    const rad = RAD[r.id]*view.k + 4;
    if ((mx-sx)**2 + (my-sy)**2 <= rad*rad) return r;
  }
  return null;
}
cv.addEventListener('wheel', e => {
  e.preventDefault();
  const rect = cv.getBoundingClientRect();
  const mx = e.clientX-rect.left, my = e.clientY-rect.top;
  const f = e.deltaY<0 ? 1.15 : 1/1.15;
  view.tx = mx - (mx-view.tx)*f; view.ty = my - (my-view.ty)*f;
  view.k *= f; drawMap();
}, {passive:false});
cv.addEventListener('mousedown', e => {
  if (e.target.closest('.fpanel')) return;
  drag = [e.clientX,e.clientY]; dragMoved = 0; cv.style.cursor='grabbing';
});
window.addEventListener('mousemove', e => {
  if (drag){
    const dx = e.clientX-drag[0], dy = e.clientY-drag[1];
    if (dragMoved > 4 || Math.abs(dx)+Math.abs(dy) > 4) dragMoved += Math.abs(dx)+Math.abs(dy);
    view.tx += dx; view.ty += dy; drag = [e.clientX,e.clientY]; drawMap(); return;
  }
  if (curTab!=='map') return;
  const rect = cv.getBoundingClientRect();
  const mx = e.clientX-rect.left, my = e.clientY-rect.top;
  if (mx<0||my<0||mx>cv.width||my>cv.height) return;
  hover = hitTest(mx, my);
  const tip = $('tip');
  if (hover){
    tip.style.display='block';
    tip.style.left = (e.clientX+14)+'px'; tip.style.top = (e.clientY+10)+'px';
    const deg = (ms.spot && degMap && degMap[hover.id] != null) ? ' · degree '+degMap[hover.id] : '';
    tip.innerHTML = '<b>'+esc(hover.t)+'</b><br>'+esc(hover.a)+' '+(hover.y||'')+' · '+hover.c+' citations'+deg+
      '<br><span class="sub">left-click: '+$('m-swap').value+' · right-click/Ctrl+click: '+($('m-swap').value==='focus'?'details':'focus')+'</span>';
  } else tip.style.display='none';
  drawMap();
});
window.addEventListener('mouseup', e => {
  if (!drag) return;
  const moved = dragMoved > 4;
  drag = null; cv.style.cursor='grab';
  if (moved || curTab!=='map') return;
  const rect = cv.getBoundingClientRect();
  const r = hitTest(e.clientX-rect.left, e.clientY-rect.top);
  if (r) act(r, e.button !== 2 && !e.ctrlKey);
});
cv.addEventListener('contextmenu', e => e.preventDefault());
cv.addEventListener('keydown', e => {
  const ns = nodesSorted();
  if (e.key === 'Escape'){
    if (ms.mode==='refs'||ms.mode==='citedby'){ if (ms.focus){ ms.focus = null; nodesSorted._c = null; resetView(); } }
    else if (ms.spot){ clearSpot(); }
    e.preventDefault(); return;
  }
  if (!ns.length) return;
  if (e.key === 'ArrowRight' || e.key === 'ArrowDown') kfocus = Math.min(ns.length-1, kfocus+1);
  else if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') kfocus = Math.max(0, kfocus-1);
  else if (e.key === 'Enter'){ if (kFocused()) act(kFocused(), !e.shiftKey); e.preventDefault(); return; }
  else return;
  e.preventDefault();
  const r = kFocused();
  if (r){
    $('map-live').textContent = r.t + ', ' + (r.y||'') + ', ' + r.c + ' citations';
    const sx = PX[r.id][0]*view.k+view.tx, sy = PX[r.id][1]*view.k+view.ty;
    if (sx < 40 || sx > cv.width-40 || sy < 40 || sy > cv.height-40) resetView();
    else drawMap();
  }
});
$('m-layout').addEventListener('change', () => {
  curLayout = $('m-layout').value;
  resetView();
});
$('m-style').addEventListener('change', drawMap);
$('m-view').addEventListener('change', () => {
  ms.mode = $('m-view').value; ms.focus = null; ms.spot = null; degMap = null; kfocus = -1;
  nodesSorted._c = null; resetView();
});
$('m-swap').addEventListener('change', drawMap);
$('m-deg').addEventListener('change', () => {
  if (ms.spot) hopFocus(ms.spot);   // recompute neighborhood at the new depth
  else drawMap();
});

// ================= DETAIL =================
function showDetail(r){
  const dd = $('detail-body');
  const nRefs = (CITES[r.id]||[]).length, nBy = (citedBy[r.id]||[]).length;
  dd.innerHTML = '<h2>'+esc(r.t)+'</h2>'+
    '<div class="sub">'+esc(r.a)+(r.s? ' — '+esc(r.s):'')+' · '+(r.y||'n.d.')+' · '+r.c+' citations</div>'+
    '<dt>DOI</dt><dd>'+(r.d? '<a href="https://doi.org/'+esc(r.d)+'" target="_blank" rel="noopener">'+esc(r.d)+'</a>' : '(none)')+'</dd>'+
    '<dt>In-corpus citations</dt><dd>'+nRefs+' references in corpus · cited by '+nBy+' corpus papers'+
    ' — <a href="#" onclick="focusInMap(\''+r.id+'\');return false;">Focus in map</a> · '+
    '<a href="#" onclick="viewRefs(\''+r.id+'\');return false;">References map</a> · '+
    '<a href="#" onclick="viewCitedBy(\''+r.id+'\');return false;">Cited-by map</a></dd>'+
    '<dt>Import source</dt><dd>'+(SRC_LABEL[r.src]||r.src||'?')+
    (r.ar? ' · <b>app-only: removed from Airtable (record-cap cut), kept here for its rules/insight</b>':'')+'</dd>'+
    '<dt>Topic cluster</dt><dd>'+(r.cl>=0? esc(clName(r.cl)) : '(no citation links)')+'</dd>'+
    '<dt>Animals</dt><dd>'+(r.an.length? r.an.map(a=>'<span class="pill">'+esc(a)+'</span>').join('') : '(untagged)')+'</dd>'+
    '<dt>Afferent types</dt><dd>'+(r.af.length? r.af.map(a=>'<span class="pill af">'+esc(a)+'</span>').join('') : '(untagged)')+'</dd>'+
    '<dt>Feedback pathways</dt><dd>'+(r.fb.length? r.fb.map(f=>'<span class="pill fb">'+esc(f)+'</span>').join('') : '(none linked)')+'</dd>'+
    '<dt>Models</dt><dd>'+esc([].concat(r.md,r.mr).filter(Boolean).join('; ')||'(none)')+'</dd>'+
    '<dt>Review coverage</dt><dd>'+esc(r.rv.join('; ')||'(none)')+'</dd>'+
    '<dt>PDF in Airtable</dt><dd>'+(r.pdf? 'yes':'no')+(r.pr? ' · prune status: '+esc(r.pr):'')+'</dd>'+
    '<dt>Curation note</dt><dd>'+(r.n? esc(r.nt) : '<i>(not yet curated)</i>')+'</dd>'+
    (r.ast? '<dt>Animal-study potential</dt><dd>'+esc(r.ast)+'</dd>':'')+
    (r.rs? '<dt>Robot / sim translation</dt><dd>'+esc(r.rs)+'</dd>':'');
  $('detail').style.display='block';
}
function hideDetail(){ $('detail').style.display='none'; }
function focusInMap(id){ hideDetail(); show('map'); hopFocus(id); }
function viewRefs(id){ ms.mode='refs'; ms.focus=id; ms.spot=null; degMap=null; $('m-view').value='refs'; hideDetail(); show('map'); }
function viewCitedBy(id){ ms.mode='citedby'; ms.focus=id; ms.spot=null; degMap=null; $('m-view').value='citedby'; hideDetail(); show('map'); }

// ================= SEARCH (online) =================
const CORPUS_DOIS = new Set(DATA.map(r => (r.d||'').toLowerCase()).filter(Boolean));
async function doSearch(){
  const q = $('s-q').value.trim();
  const res = $('s-results');
  if (!q){ res.innerHTML = '<div class="sub">Type a query first.</div>'; return; }
  res.innerHTML = '<div class="sub">Querying OpenAlex…</div>';
  const url = 'https://api.openalex.org/works?search=' + encodeURIComponent(q) +
    '&per-page=' + $('s-n').value + '&mailto=benjamin.bolen@pdx.edu' +
    '&select=id,doi,title,publication_year,authorships,cited_by_count,primary_location,abstract_inverted_index';
  let data = null;
  try {
    const r = await fetch(url);
    if (!r.ok) throw new Error('HTTP ' + r.status);
    data = await r.json();
  } catch (err){
    res.innerHTML = '<div class="warn" style="display:block">Online fetch failed ('+esc(String(err))+
      '). The OpenAlex endpoint needs internet access; everything else in this app works offline.</div>';
    return;
  }
  $('s-note').textContent = data.meta ? (data.meta.count + ' matches; showing ' + data.results.length) : '';
  if (!data.results.length){ res.innerHTML = '<div class="sub">No matches.</div>'; return; }
  res.innerHTML = data.results.map(w => {
    const doi = (w.doi||'').replace('https://doi.org/','');
    const auth = (w.authorships||[]).slice(0,4).map(a => (a.author||{}).display_name).filter(Boolean).join(', ');
    const inst = (((w.primary_location||{}).source||{})||{}).display_name || '';
    const inCorpus = doi && CORPUS_DOIS.has(doi.toLowerCase());
    let ab = '';
    if (w.abstract_inverted_index){
      const pos = {};
      for (const [word, idxs] of Object.entries(w.abstract_inverted_index)) idxs.forEach(i => pos[i]=word);
      const keys = Object.keys(pos).map(Number).sort((a,b)=>a-b);
      ab = keys.map(k => pos[k]).join(' ');
    }
    const sq = encodeURIComponent(((w.title||'') + ' ' + q).slice(0,180));
    return '<div class="srch-r"><div class="t">'+esc(w.title||'(untitled)')+
      (inCorpus? ' <span class="incorpus">already in corpus</span>':'')+'</div>'+
      '<div class="m">'+esc(auth)+' · '+(w.publication_year||'')+(inst? ' · '+esc(inst):'')+
      ' · cited by '+(w.cited_by_count||0)+(doi? ' · doi: '+esc(doi):'')+'</div>'+
      (ab? '<div class="m">'+esc(ab.slice(0,320))+(ab.length>320?'…':'')+'</div>':'')+
      '<div style="margin-top:4px">'+
      (doi? '<a href="https://doi.org/'+esc(doi)+'" target="_blank" rel="noopener">doi.org</a>':'')+
      '<a href="https://scholar.google.com/scholar?q='+sq+'" target="_blank" rel="noopener">Google Scholar ↗</a>'+
      '<a href="https://pubmed.ncbi.nlm.nih.gov/?term='+sq+'" target="_blank" rel="noopener">PubMed ↗</a>'+
      '<a href="#" onclick="openWoS(\''+esc(String(q).replace(/'/g,""))+'\');return false;">Web of Science ↗</a>'+
      '<a href="'+esc(w.id)+'" target="_blank" rel="noopener">OpenAlex ↗</a>'+
      '<a href="#" onclick="suggestImport(\''+esc((doi||w.id).replace(/'/g,""))+'\');return false;">⤓ suggest for corpus</a>'+
      '</div></div>';
  }).join('');
}
function openWoS(q){
  try { navigator.clipboard.writeText(q); } catch (e) {}
  window.open('https://www.webofscience.com/wos/woscc/basic-search', '_blank', 'noopener');
}
function suggestImport(key){
  const blob = new Blob([JSON.stringify({suggest: key, queued_by: 'SADb Explorer', 
    date: new Date().toISOString()}, null, 1)], {type: 'application/json'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'sadb_suggest_' + key.replace(/[^a-z0-9]/gi,'_').slice(0,40) + '.json';
  a.click(); URL.revokeObjectURL(a.href);
}
$('s-q').addEventListener('keydown', e => { if (e.key === 'Enter') doSearch(); });

// ================= TABS =================
let curTab = 'table';
function show(t){
  curTab = t;
  ['table','pivot','map','search'].forEach(x => {
    $('view-'+x).style.display = x===t ? '' : 'none';
    $('bar-'+x).style.display = x===t ? 'flex' : 'none';
    $('tab-'+x).className = x===t ? 'on' : '';
  });
  if (t==='table') renderTable();
  if (t==='pivot') renderPivot();
  if (t==='map'){ resetView(); }
}

// ================= INIT =================
(function init(){
  const opts = (sel, vals) => vals.forEach(v => sel.add(new Option(v, v)));
  const uniq = f => [...new Set(DATA.flatMap(f))].filter(Boolean).sort();
  opts($('f-an'), uniq(r => r.an));
  opts($('f-fb'), uniq(r => r.fb));
  opts($('f-af'), uniq(r => r.af));
  opts($('f-src'), uniq(r => [SRC_LABEL[r.src]||r.src]));
  $('f-src').addEventListener('change', () => {
    const lbl = $('f-src').value;
    state.src = Object.keys(SRC_LABEL).find(k => (SRC_LABEL[k]||k)===lbl) || lbl || '';
    renderTable();
  });
  Object.entries(DIMS).forEach(([k,d]) => {
    $('p-row').add(new Option(d.label, k));
    $('p-col').add(new Option(d.label, k));
  });
  $('p-row').value = 'cluster'; $('p-col').value = 'decade';
  ['p-row','p-col','p-sr','p-sc'].forEach(id => $(id).addEventListener('change', renderPivot));
  const withNotes = DATA.filter(r=>r.n).length, withPdf = DATA.filter(r=>r.pdf).length;
  const appOnly = DATA.filter(r=>r.ar).length;
  const ys = DATA.map(r=>r.y||9999).filter(v=>v<9999);
  $('stats').textContent = ' — '+DATA.length+' papers ('+appOnly+' app-only, removed from Airtable) · '+
    withNotes+' curated · '+withPdf+' with PDF · '+
    Math.min(...ys)+'–'+Math.max(...DATA.map(r=>r.y||0));
  renderTable();
})();
</script>
</body>
</html>
"""

html = (HTML.replace("__DATA__", payload)
            .replace("__CITES__", cites_payload)
            .replace("__CLUSTER__", labels_payload))
out = os.path.join(HERE, "sadb_app.html")
with open(out, "w", encoding="utf-8") as fh:
    fh.write(html)
print("wrote", out, f"({len(html)/1024:.0f} KB, {len(data)} papers, "
      f"{sum(len(v) for v in cites.values())} directed citation pairs, "
      f"{len(cl_labels)} named clusters) — built {datetime.date.today().isoformat()}")
