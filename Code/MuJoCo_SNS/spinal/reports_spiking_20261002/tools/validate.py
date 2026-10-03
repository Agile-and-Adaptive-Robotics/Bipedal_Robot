"""Validate an .asim after surgery + structural-diff it against its original.

Checks:
  1. whole-file XML parse
  2. every <DiagramXml><![CDATA[ ... ]]> page parses as XML on its own
  3. structural sets unchanged vs the original: neuron IDs, connexion tuples,
     synapse-type IDs, chart DataColumn TargetIDs, stimulus targets
  4. converted expectations: all neuron thresholds < 0; no NonSpikingChemical types;
     every connexion's SynapseTypeID exists
  5. SimEndTime > every chart EndTime (flush trap)
Usage: python validate.py <original.asim> <converted.asim>
"""
import io
import sys
import re
import xml.etree.ElementTree as ET

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def txt(e, t):
    x = e.find(t)
    return x.text if x is not None else None


def load(path):
    try:
        root = ET.parse(path).getroot()
        print(f"  [OK] whole-file XML parse: {path}")
        return root
    except ET.ParseError as e:
        print(f"  [FAIL] whole-file parse: {e}")
        return None


def cdata_pages(path):
    s = open(path, encoding="utf-8").read()
    pages = re.findall(r"<DiagramXml><!\[CDATA\[(.*?)\]\]></DiagramXml>", s, flags=re.S)
    ok = 0
    for i, p in enumerate(pages):
        try:
            ET.fromstring("<root>" + p + "</root>")
            ok += 1
        except ET.ParseError as e:
            print(f"  [FAIL] page CDATA #{i}: {e}")
    print(f"  [{'OK' if ok == len(pages) else 'FAIL'}] page CDATAs parsed: {ok}/{len(pages)}")
    return ok == len(pages)


def structure(root):
    neurons = {txt(n, "ID"): (txt(n, "Name"), txt(n, "InitialThresh"),
                              txt(n, "RelativeAccom")) for n in root.iter("Neuron")}
    conn = [(txt(c, "SourceID"), txt(c, "TargetID"), txt(c, "SynapseTypeID"))
            for c in root.iter("Connexion")]
    conn_kind = {txt(c, "SynapseTypeID"): txt(c, "Type") for c in root.iter("Connexion")}
    # section membership of each type id
    sect_of = {}
    for sec, kind in (("SpikingSynapses", "0"), ("NonSpikingSynapses", "1"),
                      ("ElectricalSynapses", "2")):
        for e in root.iter(sec):
            for t in e.iter("SynapseType"):
                sect_of[txt(t, "ID")] = sec
    types = {txt(t, "ID"): txt(t, "Type") for t in root.iter("SynapseType")}
    chartcols = [(txt(ch, "Name"), txt(c, "ColumnName"), txt(c, "TargetID"))
                 for ch in root.iter("DataChart") for c in ch.iter("DataColumn")]
    stims = [(txt(s, "Name"), txt(s, "TargetNodeID")) for e in root.iter("ExternalStimuli") for s in e]
    simend = [e.text for e in root.iter("SimEndTime")]
    charts = [(txt(ch, "Name"), txt(ch, "EndTime")) for ch in root.iter("DataChart")]
    return dict(neurons=neurons, conn=conn, types=types, chartcols=chartcols,
                stims=stims, simend=simend, charts=charts,
                conn_kind=conn_kind, sect_of=sect_of)


def main():
    orig, conv = sys.argv[1], sys.argv[2]
    print(f"== validate {conv} vs {orig}")
    ro = load(orig)
    rc = load(conv)
    if ro is None or rc is None:
        sys.exit(1)
    if not cdata_pages(conv):
        sys.exit(1)
    so, sc = structure(ro), structure(rc)
    ok = True
    if set(so["neurons"]) != set(sc["neurons"]):
        print("  [FAIL] neuron ID sets differ")
        ok = False
    else:
        print(f"  [OK] neuron ID sets identical ({len(sc['neurons'])})")
    if so["conn"] != sc["conn"]:
        print("  [FAIL] connexion tuples differ")
        ok = False
    else:
        print(f"  [OK] connexions identical ({len(sc['conn'])})")
    if set(so["types"]) != set(sc["types"]):
        print("  [FAIL] synapse-type ID sets differ")
        ok = False
    else:
        print(f"  [OK] synapse-type ID sets identical ({len(sc['types'])})")
    if so["chartcols"] != sc["chartcols"]:
        print("  [FAIL] chart columns differ")
        ok = False
    else:
        print(f"  [OK] chart columns identical ({len(sc['chartcols'])})")
    if so["stims"] != sc["stims"]:
        print("  [FAIL] stimuli differ")
        ok = False
    else:
        print(f"  [OK] stimuli identical ({len(sc['stims'])})")
    # converted expectations
    nonneg = [v[0] for v in sc["neurons"].values() if float(v[1]) >= 0]
    if nonneg:
        print(f"  [FAIL] neurons still non-spiking (thr>=0): {nonneg[:5]}...")
        ok = False
    else:
        print("  [OK] all neuron thresholds < 0 (spiking regime)")
    ns_types = [i for i, t in sc["types"].items() if t == "NonSpikingChemical"]
    if ns_types:
        print(f"  [FAIL] NonSpikingChemical types remain: {len(ns_types)}")
        ok = False
    else:
        print("  [OK] no NonSpikingChemical types remain")
    missing = {c[2] for c in sc["conn"]} - set(sc["types"])
    if missing:
        print(f"  [FAIL] connexions pointing at missing types: {len(missing)}")
        ok = False
    else:
        print("  [OK] every connexion's SynapseTypeID resolves")
    # connexion kind enum vs type section/class consistency
    bad_kind = []
    for sid, kind in sc["conn_kind"].items():
        sec = sc["sect_of"].get(sid, "?")
        expect = {"SpikingSynapses": "0", "NonSpikingSynapses": "1",
                  "ElectricalSynapses": "2"}.get(sec)
        if expect is not None and kind != expect:
            bad_kind.append((sid[:8], kind, sec))
    if bad_kind:
        print(f"  [FAIL] connexion <Type> enum vs synapse class mismatch: {bad_kind[:3]}")
        ok = False
    else:
        print("  [OK] every connexion's <Type> enum matches its synapse class")
    # chart-guid validity: every TargetID must exist among known IDs
    known = set(sc["neurons"])
    for b in rc.iter("Body"):
        if txt(b, "ID"):
            known.add(txt(b, "ID"))
    for j in rc.iter("Joint"):
        if txt(j, "ID"):
            known.add(txt(j, "ID"))
    for m in rc.iter("Muscle"):
        if txt(m, "ID"):
            known.add(txt(m, "ID"))
    badcol = [c for c in sc["chartcols"] if c[2] not in known]
    if badcol:
        print(f"  [WARN] chart columns with unknown TargetID (chart silently writes 0 bytes): {badcol[:3]}")
    else:
        print("  [OK] all chart-column TargetIDs resolve to body/neuron GUIDs")
    se = float(sc["simend"][0]) if sc["simend"] else None
    badend = [c for c in sc["charts"] if float(c[1] or 0) >= se]
    if badend:
        print(f"  [FAIL] chart EndTime >= SimEndTime: {badend}")
        ok = False
    else:
        print(f"  [OK] SimEndTime {se} exceeds all chart EndTimes")
    print("  VALIDATION " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
