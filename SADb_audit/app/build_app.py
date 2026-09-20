"""Build the single-file SADb explorer app: app/sadb_app.html.

Reads export/sadb_export.json + export/sadb_cites.json and embeds everything
into ONE self-contained HTML file (no CDN, no network — double-click to open).
Tabs: Table (search/sort/filter), Pivot (count matrix + drill-through), and a
bubble map on PUBLICATION YEAR (x) vs CITATION COUNT (y, log scale).

Map interactions (accessible + Mac friendly):
  View dropdown: All papers / References / Cited by / Research review /
  Animal studies / Models. In References/Cited-by, left-clicking a bubble makes
  that paper the focus of the map; in other modes it spotlights the paper with
  its corpus citation neighborhood. The "Left click action" selector swaps the
  primary/secondary actions. Secondary action = right-click, Ctrl+click
  (Mac trackpad), or Shift+Enter. Keyboard: Tab to the map, arrow keys move
  between bubbles, Enter = primary, Shift+Enter = secondary, Esc = back.
Bubble colors: Ben's accessible 7-color set from Code\\Matlab\\Colors.m.

Re-run after any curation batch:  myo python export_corpus.py && vos_build2.py && app/build_app.py
"""
import datetime, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
records = json.load(open(os.path.join(SAD, "export", "sadb_export.json"), encoding="utf-8"))
cites = json.load(open(os.path.join(SAD, "export", "sadb_cites.json"), encoding="utf-8"))

# --- source classification (same rule as vos_build2.py; export json lacks src) ---
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
        "md": r["models2"], "mr": r["models_ref"], "rv": r["reviews"],
        "pdf": 1 if r["has_pdf"] else 0, "n": 1 if r["has_notes"] else 0,
        "nt": r["notes"], "src": _classify(r), "c": 0, "cl": -1,
    })

layout_path = os.path.join(SAD, "export", "sadb_layout.json")
if os.path.exists(layout_path):
    _lay = json.load(open(layout_path, encoding="utf-8"))
    for r in data:
        L = _lay.get(r["id"], {})
        r["c"] = L.get("c", 0)
        r["cl"] = L.get("cl", -1)
else:
    print("WARNING: sadb_layout.json missing — citation counts will be 0")

payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
cites_payload = json.dumps(cites, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

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
 #count { color:var(--mut); font-size:12px; margin:6px 0; }
 .pivot td, .pivot th { text-align:right; }
 .pivot td.rh, .pivot th.rh { text-align:left; }
 .pivot td.cell { cursor:pointer; }
 .pivot td.cell:hover { background:#ffe9b8; }
 .pivot td.tot { font-weight:600; background:#f4f4f4; }
 #crumb { display:flex; align-items:center; gap:10px; padding:4px 0; font-size:13px; }
 #crumb button { border:1px solid var(--acc); background:#fff; color:var(--acc); padding:2px 10px; cursor:pointer; }
 #canvas { border:1px solid var(--line); background:#fff; cursor:grab; display:block; outline-offset:-3px; }
 #canvas:focus { outline:3px solid var(--acc); }
 #tip { position:fixed; display:none; background:#fff; border:1px solid var(--acc); padding:6px 9px;
        font-size:12px; max-width:340px; pointer-events:none; box-shadow:2px 2px 8px rgba(0,0,0,.15); z-index:50; }
 #detail { position:fixed; right:0; top:0; bottom:0; width:430px; background:#fff; border-left:2px solid var(--acc);
           padding:14px 16px; overflow-y:auto; display:none; z-index:40; box-shadow:-4px 0 16px rgba(0,0,0,.12); }
 #detail h2 { font-size:15px; margin:0 0 6px; } #detail .close { float:right; cursor:pointer; border:1px solid var(--line); padding:1px 8px; }
 #detail dt { font-weight:600; font-size:12px; color:var(--mut); margin-top:8px; }
 #detail dd { margin:2px 0 0; }
 .legend { display:flex; flex-wrap:wrap; gap:4px 12px; font-size:11px; color:var(--mut); margin-top:6px; }
 .legend i { display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:3px; }
 .warn { background:#fff4d5; border:1px solid #e8d49a; padding:4px 10px; font-size:12px; display:none; }
 .sr { position:absolute; width:1px; height:1px; overflow:hidden; clip:rect(0 0 0 0); }
</style>
</head>
<body>
<header>
  <h1>SADb Explorer — Sensory Afferent Database <span class="sub" id="stats"></span></h1>
  <div class="sub">Offline snapshot of the Airtable "Sensory Feedback" corpus. Rebuild: <code>myo python SADb_audit/export_corpus.py &amp;&amp; SADb_audit/vos_build2.py &amp;&amp; SADb_audit/app/build_app.py</code></div>
  <div class="tabs">
    <button id="tab-table" class="on" onclick="show('table')">Table</button>
    <button id="tab-pivot" onclick="show('pivot')">Pivot</button>
    <button id="tab-map" onclick="show('map')">Bubble map</button>
  </div>
</header>

<div class="bar" id="bar-table">
  <input type="text" id="q" placeholder="Search title / author / year / DOI…" size="38">
  <label>Animal <select id="f-an"><option value="">(any)</option></select></label>
  <label>Pathway <select id="f-fb"><option value="">(any)</option></select></label>
  <label>Source <select id="f-src"><option value="">(any)</option></select></label>
  <label>Year <input type="number" id="f-y0" placeholder="from"> – <input type="number" id="f-y1" placeholder="to"></label>
  <label><input type="checkbox" id="f-notes"> has notes</label>
  <label><input type="checkbox" id="f-pdf"> has PDF</label>
  <span class="warn" id="cellwarn"></span>
  <button onclick="clearCell()">Clear drill-through</button>
</div>

<div class="bar" id="bar-pivot" style="display:none">
  <label>Rows <select id="p-row"></select></label>
  <label>Columns <select id="p-col"></select></label>
  <label class="sub">Count = papers. Multi-valued fields (animal, pathway) count the paper once per value. Click a cell to drill through to the Table tab.</label>
</div>

<div class="bar" id="bar-map" style="display:none">
  <label>View <select id="m-view">
    <option value="all">All papers</option>
    <option value="refs">References</option>
    <option value="citedby">Cited by</option>
    <option value="review">Research review</option>
    <option value="animal">Animal studies</option>
    <option value="models">Models</option>
  </select></label>
  <label>Left click action <select id="m-swap">
    <option value="focus">Focus paper</option>
    <option value="details">Show details</option>
  </select></label>
  <label class="sub">Left-click = focus (in References/Cited-by: make that paper the map; otherwise spotlight it). Right-click, Ctrl+click (Mac) or Shift+Enter = details. Tab/arrow keys + Enter also work; Esc = back. Wheel = zoom, drag = pan.</label>
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
          aria-label="Bubble map: publication year (horizontal) versus citation count (vertical). Use arrow keys to move between bubbles, Enter for the primary action, Shift+Enter for details, Escape to go back."></canvas>
  <div class="sr" id="map-live" aria-live="polite"></div>
  <div class="legend" id="maplegend"></div>
</main>

<div id="tip"></div>
<div id="detail"><span class="close" onclick="hideDetail()">✕</span><div id="detail-body"></div></div>

<script type="application/json" id="sadb-data">__DATA__</script>
<script type="application/json" id="sadb-cites">__CITES__</script>
<script>
const DATA = JSON.parse(document.getElementById('sadb-data').textContent);
const CITES = JSON.parse(document.getElementById('sadb-cites').textContent);
/* Ben's accessible palette from Code\Matlab\Colors.m */
const PAL = ['#FFD700','#FFB14E','#FA8775','#EA5F94','#CD34B5','#9D02D7','#0000FF'];
const GREY = '#c9c9c9';
const SRC_LABEL = {originals99:'original 99', demo50:'demo 50', rest383:'rest import 383',
                   digest20:'RR digest 20', task5import:'task5 import'};
const byId = {}; DATA.forEach(r => byId[r.id] = r);
const citedBy = {};
for (const [p, refs] of Object.entries(CITES)) refs.forEach(t => { (citedBy[t] = citedBy[t] || []).push(p); });
const ANIMALS = [...new Set(DATA.flatMap(r => r.an))].sort();

const state = { q:'', an:'', fb:'', src:'', y0:null, y1:null, notes:false, pdf:false,
                sort:{k:'c', dir:-1}, cell:null };
const $ = id => document.getElementById(id);
function esc(s){ return (s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;'); }

// ================= TABLE =================
function hay(r){ return (r.t+' '+r.a+' '+r.s+' '+r.y+' '+r.d).toLowerCase(); }
function matches(r){
  const st = state;
  if (st.q && !hay(r).includes(st.q.toLowerCase())) return false;
  if (st.an && !r.an.includes(st.an)) return false;
  if (st.fb && !r.fb.includes(st.fb)) return false;
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
  source:  { label:'Import source', vals:r=>[SRC_LABEL[r.src]||r.src||'?'] },
  decade:  { label:'Decade', vals:r=>[r.y? Math.floor(r.y/10)*10+'s' : 'no year'] },
  author:  { label:'Primary author', vals:r=>[r.a||'(none)'] },
  animal:  { label:'Animal', multi:true, vals:r=> r.an.length? r.an : ['(untagged)'] },
  pathway: { label:'Feedback pathway', multi:true, vals:r=> r.fb.length? r.fb : ['(no pathway link)'] },
  cluster: { label:'Topic cluster', vals:r=>[r.cl>=0? 'Cluster '+(r.cl+1) : 'no citation links'] },
  notes:   { label:'Notes', vals:r=>[r.n? 'has notes':'no notes'] },
  pdf:     { label:'PDF', vals:r=>[r.pdf? 'has PDF':'no PDF'] },
};
function dimHas(r, key, val){ return DIMS[key].vals(r).includes(val); }

function renderTable(){
  const rows = filtered();
  $('count').textContent = rows.length + ' of ' + DATA.length + ' papers' +
    (state.cell ? '  [drill-through active]' : '');
  const tb = $('tbody'); tb.innerHTML = '';
  const frag = document.createDocumentFragment();
  rows.slice(0, 400).forEach(r => {
    const tr = document.createElement('tr'); tr.className = 'row';
    const tags = [];
    if (r.n) tags.push('<span class="pill yes">notes</span>');
    else tags.push('<span class="pill no">no notes</span>');
    if (r.pdf) tags.push('<span class="pill yes">PDF</span>');
    r.an.slice(0,4).forEach(a => tags.push('<span class="pill">'+a+'</span>'));
    r.fb.slice(0,3).forEach(f => tags.push('<span class="pill fb">'+f+'</span>'));
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
['q','f-an','f-fb','f-src','f-y0','f-y1','f-notes','f-pdf'].forEach(id => {
  $(id).addEventListener(id==='q' ? 'input' : 'change', () => {
    state.q = $('q').value; state.an = $('f-an').value; state.fb = $('f-fb').value;
    state.src = $('f-src').value;
    state.y0 = $('f-y0').value ? +$('f-y0').value : null;
    state.y1 = $('f-y1').value ? +$('f-y1').value : null;
    state.notes = $('f-notes').checked; state.pdf = $('f-pdf').checked;
    renderTable(); if (curTab==='map') drawMap();
  });
});
function clearCell(){ state.cell = null; $('cellwarn').style.display='none'; renderTable(); }

// ================= PIVOT =================
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
  const rl = [...rVals.entries()].sort((a,b)=>b[1]-a[1]).slice(0,40);
  const cl = [...cVals.entries()].sort((a,b)=>b[1]-a[1]).slice(0,24);
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
  if (rVals.size>40 || cVals.size>24) h += '<p class="sub">Showing top 40 row / 24 column values by total.</p>';
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
const POS = {}, RAD = {};
DATA.forEach(r => {
  const jh = hashId(r.id);
  POS[r.id] = [
    (r.y || 1970) + ((jh % 1000)/1000 - 0.5) * 1.8,
    Math.log10((r.c||0) + 1) + (((jh>>>10) % 1000)/1000 - 0.5) * 0.16,
  ];
  RAD[r.id] = Math.min(24, 1.7 + 5.5*Math.pow((r.c||0)/MAXC, 0.25)*Math.pow(MAXC, 0.12));
});
function bx(x){ return PAD.l + (x - X0)/(X1 - X0) * (cv.width - PAD.l - PAD.r); }
function by(y){ return cv.height - PAD.b - (y - Y0)/(Y1 - Y0) * (cv.height - PAD.t - PAD.b); }
function radOf(r){ return RAD[r.id]; }

function nodesNow(){
  if ((ms.mode==='refs' || ms.mode==='citedby') && ms.focus){
    const nb = ms.mode==='refs' ? (CITES[ms.focus]||[]) : (citedBy[ms.focus]||[]);
    const set = new Set([ms.focus, ...nb]);
    return DATA.filter(r => set.has(r.id));
  }
  return DATA;
}
function colorOf(r){
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
  // axes (fixed frame)
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

  // bubbles: large first so small ones sit on top and stay clickable
  const shown = nodesNow();
  const searchSet = state.q ? new Set(filtered().map(r=>r.id)) : null;
  const inSpot = ms.spot ? new Set([ms.spot, ...(CITES[ms.spot]||[]), ...(citedBy[ms.spot]||[])]) : null;
  const ordered = [...shown].sort((a,b) => RAD[b.id]-RAD[a.id]);
  ctx.setTransform(view.k,0,0,view.k,view.tx,view.ty);
  ordered.forEach(r => {
    const [X,Y] = POS[r.id]; const sx = bx(X), sy = by(Y);
    ctx.beginPath(); ctx.arc(sx, sy, RAD[r.id], 0, 6.2832);
    let alpha = 0.85;
    if (searchSet && !searchSet.has(r.id)) alpha = 0.10;
    if (inSpot && !inSpot.has(r.id)) alpha = 0.10;
    ctx.globalAlpha = alpha;
    ctx.fillStyle = colorOf(r); ctx.fill();
    if (r.id === (ms.focus||'')){ ctx.globalAlpha = 1; ctx.lineWidth = 2.5/view.k; ctx.strokeStyle = '#111'; ctx.stroke(); }
    if (r.id === (ms.spot||'')){ ctx.globalAlpha = 1; ctx.lineWidth = 2.5/view.k; ctx.strokeStyle = '#111'; ctx.stroke(); }
    if (r === hover || (kfocus >= 0 && r === kFocused())) {
      ctx.globalAlpha = 1; ctx.lineWidth = 2/view.k; ctx.strokeStyle = '#000'; ctx.setLineDash([4/view.k,3/view.k]); ctx.stroke(); ctx.setLineDash([]);
    }
  });
  ctx.globalAlpha = 1;
  ctx.setTransform(1,0,0,1,0,0);
  renderLegend(); updateCrumb();
}
function kFocused(){ return kfocus >= 0 ? nodesSorted()[kfocus] : null; }
function nodesSorted(){ // keyboard order: by year then citations
  const ns = nodesNow();
  if (!nodesSorted._c || nodesSorted._c !== ns){
    nodesSorted._list = [...ns].sort((a,b)=> (a.y||0)-(b.y||0) || a.c-b.c);
    nodesSorted._c = ns;
  }
  return nodesSorted._list;
}

function renderLegend(){
  const m = ms.mode; let items = [];
  if (m==='review') items = [['has Review Papers link', PAL[5]], ['not tagged as review', GREY]];
  else if (m==='models') items = [['is a model (Models 2)', PAL[4]], ['references models', PAL[1]], ['neither', GREY]];
  else if (m==='animal'){
    const used = {}; DATA.forEach(r => { if (r.an.length) used[r.an[0]] = (used[r.an[0]]||0)+1; });
    items = Object.entries(used).sort((a,b)=>b[1]-a[1]).slice(0,14)
      .map(([a,n]) => [a+' ('+n+')', PAL[ANIMALS.indexOf(a)%PAL.length]]);
  } else {
    const used = {}; DATA.forEach(r => { if (r.cl>=0) used[r.cl] = (used[r.cl]||0)+1; });
    items = Object.keys(used).sort((a,b)=>used[b]-used[a]).slice(0,14)
      .map(c => ['Cluster '+(+c+1)+' ('+used[c]+')', PAL[c%PAL.length]]);
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
      esc(r.t.slice(0,80));
  } else {
    const lbl = {all:'All papers', review:'Research review', animal:'Animal studies', models:'Models'}[ms.mode];
    c.innerHTML = '<span class="sub">'+lbl+' — click a bubble'+
      (ms.mode==='all' ? '; switch View to References or Cited by, then click a paper to see its citation neighborhood.' : '.')+'</span>';
  }
  const back = $('map-back');
  if (back) back.onclick = () => {
    if (ms.mode==='refs'||ms.mode==='citedby'){ ms.focus = null; } else { ms.spot = null; }
    kfocus = -1; nodesSorted._c = null; drawMap();
  };
}
function resetView(){
  sizeCanvas();               // fit must be computed in the live canvas frame
  view = {k:1, tx:0, ty:0};
  const ns = nodesNow();
  const xs = ns.map(r=>bx(POS[r.id][0])), ys = ns.map(r=>by(POS[r.id][1]));
  const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  view.k = Math.max(0.4, Math.min(cv.width/(x1-x0+80), cv.height/(y1-y0+80), 4));
  view.tx = cv.width/2 - (x0+x1)/2*view.k; view.ty = cv.height/2 - (y0+y1)/2*view.k;
  drawMap();
}
function primaryAct(r){
  if (ms.mode==='refs'){ ms.focus = r.id; ms.spot = null; kfocus = -1; nodesSorted._c = null; resetView(); return; }
  if (ms.mode==='citedby'){ ms.focus = r.id; ms.spot = null; kfocus = -1; nodesSorted._c = null; resetView(); return; }
  ms.spot = r.id; nodesSorted._c = null; drawMap();
}
function secondaryAct(r){ showDetail(r); }
function act(r, primary){
  const doFocus = $('m-swap').value==='focus' ? primary : !primary;
  if (doFocus) primaryAct(r); else secondaryAct(r);
}

function hitTest(mx, my){
  const ordered = [...nodesNow()].sort((a,b)=>RAD[a.id]-RAD[b.id]); // small first = topmost first
  for (const r of ordered){
    const sx = bx(POS[r.id][0])*view.k + view.tx, sy = by(POS[r.id][1])*view.k + view.ty;
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
cv.addEventListener('mousedown', e => { drag = [e.clientX,e.clientY]; dragMoved = 0; cv.style.cursor='grabbing'; });
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
    tip.innerHTML = '<b>'+esc(hover.t)+'</b><br>'+esc(hover.a)+' '+(hover.y||'')+' · '+hover.c+' citations'+
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
    else if (ms.spot){ ms.spot = null; nodesSorted._c = null; drawMap(); }
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
    const [X,Y] = POS[r.id]; const sx = bx(X)*view.k+view.tx, sy = by(Y)*view.k+view.ty;
    if (sx < 40 || sx > cv.width-40 || sy < 40 || sy > cv.height-40) resetView();
    else drawMap();
  }
});
$('m-view').addEventListener('change', () => {
  ms.mode = $('m-view').value; ms.focus = null; ms.spot = null; kfocus = -1;
  nodesSorted._c = null; resetView();
});
$('m-swap').addEventListener('change', drawMap);

// ================= DETAIL =================
function showDetail(r){
  const dd = $('detail-body');
  const nRefs = (CITES[r.id]||[]).length, nBy = (citedBy[r.id]||[]).length;
  dd.innerHTML = '<h2>'+esc(r.t)+'</h2>'+
    '<div class="sub">'+esc(r.a)+(r.s? ' — '+esc(r.s):'')+' · '+(r.y||'n.d.')+' · '+r.c+' citations</div>'+
    '<dt>DOI</dt><dd>'+(r.d? '<a href="https://doi.org/'+esc(r.d)+'" target="_blank" rel="noopener">'+esc(r.d)+'</a>' : '(none)')+'</dd>'+
    '<dt>In-corpus citations</dt><dd>'+nRefs+' references in corpus · cited by '+nBy+' corpus papers'+
      ' — <a href="#" onclick="viewRefs(\''+r.id+'\');return false;">References map</a> · '+
      '<a href="#" onclick="viewCitedBy(\''+r.id+'\');return false;">Cited-by map</a></dd>'+
    '<dt>Import source</dt><dd>'+(SRC_LABEL[r.src]||r.src||'?')+'</dd>'+
    '<dt>Animals</dt><dd>'+(r.an.length? r.an.map(a=>'<span class="pill">'+esc(a)+'</span>').join('') : '(untagged)')+'</dd>'+
    '<dt>Feedback pathways</dt><dd>'+(r.fb.length? r.fb.map(f=>'<span class="pill fb">'+esc(f)+'</span>').join('') : '(none linked)')+'</dd>'+
    '<dt>Models</dt><dd>'+esc([].concat(r.md,r.mr).filter(Boolean).join('; ')||'(none)')+'</dd>'+
    '<dt>Review coverage</dt><dd>'+esc(r.rv.join('; ')||'(none)')+'</dd>'+
    '<dt>PDF in Airtable</dt><dd>'+(r.pdf? 'yes':'no')+'</dd>'+
    '<dt>Curation note</dt><dd>'+(r.n? esc(r.nt) : '<i>(not yet curated)</i>')+'</dd>';
  $('detail').style.display='block';
}
function hideDetail(){ $('detail').style.display='none'; }
function viewRefs(id){ ms.mode='refs'; ms.focus=id; ms.spot=null; $('m-view').value='refs'; hideDetail(); show('map'); }
function viewCitedBy(id){ ms.mode='citedby'; ms.focus=id; ms.spot=null; $('m-view').value='citedby'; hideDetail(); show('map'); }

// ================= TABS =================
let curTab = 'table';
function show(t){
  curTab = t;
  ['table','pivot','map'].forEach(x => {
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
  ['p-row','p-col'].forEach(id => $(id).addEventListener('change', renderPivot));
  const withNotes = DATA.filter(r=>r.n).length, withPdf = DATA.filter(r=>r.pdf).length;
  const ys = DATA.map(r=>r.y||9999).filter(v=>v<9999);
  $('stats').textContent = ' — '+DATA.length+' papers · '+withNotes+' curated · '+withPdf+' with PDF · '+
    Math.min(...ys)+'–'+Math.max(...DATA.map(r=>r.y||0));
  renderTable();
})();
</script>
</body>
</html>
"""

html = (HTML.replace("__DATA__", payload)
            .replace("__CITES__", cites_payload))
out = os.path.join(HERE, "sadb_app.html")
with open(out, "w", encoding="utf-8") as fh:
    fh.write(html)
print("wrote", out, f"({len(html)/1024:.0f} KB, {len(data)} papers, "
      f"{sum(len(v) for v in cites.values())} directed citation pairs) — built {datetime.date.today().isoformat()}")
