r"""Milestone 5: AFFERENT FEEDBACK on the M4 split-RG W2L net, wired according
to Ben's circuit-rules file (canonical copy spinal\ben_rules_20260924.json =
Neuromechanical_Models\Mujoco_SNS_models\Circuit_rules_CONNECTOME_md__connectome.json).

The M5 ask fixes the priority order and the gains (quoted verbatim):
  (1) HEEL = ipsilateral stance-phase reset AT THE PF LAYER
      (heel IN -> InE/InF + the PF-layer INs, g 0.5)
  (2) TOE  = dorsiflexion inhibition ONLY
      (toe IN -> TOEDF IN -> PF-Dorsiflexion inh, g 5)
  (3) Ib load unchanged (PF-E 0.5, RG-E gate 0.1) from muscle force sensors

EDGE-BY-EDGE TRANSCRIPTION from ben_rules_20260924.json (label -> our M4 label):

  heel (node "heel: heel", SN-heel) — 5 drawn out-edges, all exc g 0.5:
    heel -> IN-InE_2806        = ipsi  "S RG ext IN"  (inhibits ipsi RG-F:
                               prolongs ipsi stance)                    g 0.5
    heel -> IN-InF_2813        = contra "S' RG flx IN" (inhibits contra RG-E
                               [n610 "HC-RG-F_Contralateral" is type HC-RG-E]:
                               releases contra flexion)                 g 0.5
    heel -> IN-PF_2810         = ipsi  "S Hip PF ext IN"  (the E-side laminate
    heel -> IN-PF_2812           INs of the hip/knee PF layers -> biases PF
                               output toward extension = the reset AT THE
                               PF LAYER)                                g 0.5
    heel -> IN-PF_2814         = the ankle PF-E IN: our 2023 net has NO ankle
                               PF layer (ankle MNs ride the knee PFs, M3) ->
                               N/A here; the knee PF-E IN carries the reset
                               the ankle path would have received. DOCUMENTED
                               REDUCTION, not a silent drop.

  toe (node "toe: toe", SN-toe) -> "toe: toe IN" (IN-C) exc 1 ->
       "IN-PF_dorsiflexion_inhibit" exc 5 -> "HC-PF-Dorsiflexion" inh 2.749.
    Our 2023 net has no dedicated dorsiflexion PF half-center (the DF MN is
    driven by the knee PF flx through pf_to_mn), so the terminal inhibition
    lands on the DORSIFLEXION MN POOL "S Ank MN flx" — the DF output element
    of our reduced net. DOCUMENTED REDUCTION (the ask's own short form says
    "toe IN -> TOEDF IN -> PF-Dorsiflexion inh, g 5": the inhibitory target is
    the dorsiflexion drive, which in this net is the DF MN).

  Ib load (node "ib_load: Ib grp", SN-Ib) — drawn out-edges exc 0.5 onto
    ipsi RG-E ("S RG ext"), ipsi InE ("S RG ext IN"), ipsi PF-E ("S Hip PF
    ext", "S Knee PF ext"). The drawing's 4th edge "RG-E -> HC-PF-E 0.1" is
    the NORMAL RG->PF drive (present in our net as pf_drive; the M4 report
    documents the 0.1-drawing vs 2.4-stack scale difference) — not re-added.
    Ankle PF-E: N/A (no ankle PF layer), same reduction as heel.

  Driven by: heel/toe from SCRIPTED stance-phase pulses matching the stepping
  phase (the runner's HEEL_c/TOE_c contact-port pattern; real contact arrives
  at milestone 6), Ib grp from the ipsilateral EXTENSOR-muscle force sensors
  (Golgi-like: mean normalized force of the side's hip/knee/ankle extension
  actuators; the R hip crossing of M3 MUSCLE_MAP is respected — the R-side
  extensor actuator set is {hip_R_flx, knee_R_ext, ankle_R_ext}).

CONDITIONAL TOPOLOGY (lab law): each afferent family's neurons AND synapses
are built only when its switch is on; the PORT-load source edges are not
synapses (zero current through them leaves the graph unchanged).

CENSUS (asserted, families heel+toe+ib all on):
  97 neurons / 188 synapses / 12 inputs / 12 outputs
  (M4 coupled baseline 87/166/6/12; +10 sensor/interneurons, +22 synapses:
   heel_rge_in 4, heel_pf 4, toe_in 2, toe_df 2, toe_df_mn 2, ib_load 8)

USAGE
    net = build()                       # all three families on, Ben's gains
    u = net.make_inputs()               # len 12
    u[net.input_index("PORT heel L")] = 4.0   # nA into the L heel sensor
    u[net.input_index("PORT Ib R")] = current # nA into the R Ib group cell
    V = net.step(u)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import build_w2l_split_net as SPLIT  # noqa: E402  (M4 net we extend)

# additive type-table keys (types that never occur in the M4 template, so the
# M4 build is untouched): contact-encoder membrane = w2l_cpg TAU["contact"].
SPLIT._TYPE_TAU.update({"SN-heel": "contact", "SN-toe": "contact"})

DT = SPLIT.DT
E_HI = SPLIT.E_HI
MUSCLE_MAP = SPLIT.MUSCLE_MAP
W2L_GAINS = SPLIT.W2L_GAINS

# Ben's gains (the ask's literal numbers) as default gain keys.
AFF_GAINS = dict(aff_heel=0.5,        # heel -> InE/InF + PF-E INs (rules g 0.5)
                 aff_toe_in=1.0,      # toe sensor -> toe IN (rules g 1)
                 aff_toe_df=5.0,      # toe IN -> TOEDF IN (rules g 5)
                 aff_toe_df_mn=2.749, # TOEDF IN -> DF drive inh (rules 2.749)
                 aff_ib=0.5)          # Ib grp -> RG-E/InE/PF-E (rules g 0.5)

SIDES = ("L", "R")
CONTRA = {"L": "R", "R": "L"}

# ---- edge specs (from, to, sign, tag) — ben_rules_20260924.json verbatim,
# label-mapped to the M4 template. f-strings over S = side, C = contra.
HEEL_EDGES = []
for S in SIDES:
    C = CONTRA[S]
    HEEL_EDGES += [
        (f"heel {S}", f"{S} RG ext IN", "exc", "heel_rge_in"),   # ipsi InE
        (f"heel {S}", f"{C} RG flx IN", "exc", "heel_rge_in"),   # contra InF
        (f"heel {S}", f"{S} Hip PF ext IN", "exc", "heel_pf"),
        (f"heel {S}", f"{S} Knee PF ext IN", "exc", "heel_pf"),
    ]

TOE_EDGES = []
for S in SIDES:
    TOE_EDGES += [
        (f"toe {S}", f"toe IN {S}", "exc", "toe_in"),
        (f"toe IN {S}", f"TOEDF IN {S}", "exc", "toe_df"),
        (f"TOEDF IN {S}", f"{S} Ank MN flx", "inh", "toe_df_mn"),
    ]

IB_EDGES = []
for S in SIDES:
    IB_EDGES += [
        (f"Ib grp {S}", f"{S} RG ext", "exc", "ib_load"),
        (f"Ib grp {S}", f"{S} RG ext IN", "exc", "ib_load"),
        (f"Ib grp {S}", f"{S} Hip PF ext", "exc", "ib_load"),
        (f"Ib grp {S}", f"{S} Knee PF ext", "exc", "ib_load"),
    ]

AFF_NODES = []   # (label, type) — built per family below
for S in SIDES:
    AFF_NODES += [
        (f"heel {S}", "SN-heel"),
        (f"toe {S}", "SN-toe"),
        (f"toe IN {S}", "IN-C"),
        (f"TOEDF IN {S}", "IN-PF"),
        (f"Ib grp {S}", "SN-Ib"),
    ]


def make_aff_template(split_tpl: dict, heel: bool = True, toe: bool = True,
                      ib: bool = True) -> dict:
    """M4 split template -> M5 afferented template (Ben's rules edges).

    Family switches drop their neurons AND synapses entirely (conditional
    topology: zero-gain synapses still change BLAS summation order)."""
    nodes = [dict(n) for n in split_tpl["nodes"]]
    edges = [dict(e) for e in split_tpl["edges"]]
    labels = {n["label"] for n in nodes}

    add_edges: list[tuple] = []
    if heel:
        add_edges += HEEL_EDGES
    if toe:
        add_edges += TOE_EDGES
    if ib:
        add_edges += IB_EDGES
    if not add_edges:
        return {"nodes": nodes, "edges": edges}

    used = {f for f, _, _, _ in add_edges} | {t for _, t, _, _ in add_edges}
    for lab, typ in AFF_NODES:
        if lab in used:
            assert lab not in labels, f"{lab} already in template"
            nodes.append({"label": lab, "type": typ})
            labels.add(lab)
    # current source ports on the three sensor neurons per side (the runner's
    # contact-port pattern; scripted in the gate, real contact at milestone 6)
    for S in SIDES:
        for fam in ("heel", "toe", "Ib"):
            if fam == "heel" and not heel:
                continue
            if fam == "toe" and not toe:
                continue
            if fam == "Ib" and not ib:
                continue
            port = f"PORT {fam} {S}"
            tgt = f"Ib grp {S}" if fam == "Ib" else f"{fam} {S}"
            nodes.append({"label": port, "type": "PORT-load"})
            labels.add(port)
            edges.append({"from": port, "to": tgt, "sign": "exc",
                          "tag": "aff_port"})
    edges += [{"from": f, "to": t, "sign": s, "tag": g}
              for f, t, s, g in add_edges]

    # ---- Ben's-rules sanity asserts (structure is pinned, not assumed)
    for S in SIDES:
        C = CONTRA[S]
        for f, t in ((f"heel {S}", f"{S} RG ext IN"),
                     (f"heel {S}", f"{C} RG flx IN"),
                     (f"heel {S}", f"{S} Hip PF ext IN"),
                     (f"heel {S}", f"{S} Knee PF ext IN")):
            assert any(e["from"] == f and e["to"] == t and
                       e["sign"] == "exc" for e in edges), f"missing {f}->{t}"
        assert any(e["from"] == f"toe IN {S}" and e["to"] == f"TOEDF IN {S}"
                   and e["tag"] == "toe_df" for e in edges)
        assert any(e["from"] == f"TOEDF IN {S}" and e["to"] == f"{S} Ank MN flx"
                   and e["sign"] == "inh" for e in edges)
    return {"nodes": nodes, "edges": edges}


@dataclass
class W2LAffNet(SPLIT.W2LSplitNet):
    """M4 split net + Ben's-rules afferents (heel / toe / Ib load)."""
    heel_on: bool = True
    toe_on: bool = True
    ib_on: bool = True

    def __post_init__(self):
        # the parent __post_init__ builds self.G = W2L_GAINS + self.gains, so
        # inject Ben's afferent gains INTO self.gains (caller overrides win).
        self.gains = {**AFF_GAINS, **(self.gains or {})}
        super().__post_init__()

    def _gain(self, tag: str, sign: str, e: dict, n_types: dict) -> float:
        G = self.G
        if tag == "heel_rge_in" or tag == "heel_pf":
            return G["aff_heel"]
        if tag == "toe_in":
            return G["aff_toe_in"]
        if tag == "toe_df":
            return G["aff_toe_df"]
        if tag == "toe_df_mn":
            return G["aff_toe_df_mn"]
        if tag == "ib_load":
            return G["aff_ib"]
        return super()._gain(tag, sign, e, n_types)

    def _assert_complete(self, edges, n_types):
        expected = {}
        for e in edges:
            if e["tag"] == "mn_to_muscle":
                continue
            if n_types[e["from"]] == "PORT-load":
                continue
            expected[e["tag"]] = expected.get(e["tag"], 0) + 1
        # family-dependent census: M4 coupled baseline 87/166/6 + per side
        # (heel: 1 neuron + 4 syn, toe: 3 neurons + 3 syn, Ib: 1 neuron + 4 syn)
        n_syn = 166 + (8 if self.heel_on else 0) + (6 if self.toe_on else 0) \
            + (8 if self.ib_on else 0)
        n_neu = 87 + (2 if self.heel_on else 0) + (6 if self.toe_on else 0) \
            + (2 if self.ib_on else 0)
        n_in = 6 + (2 if self.heel_on else 0) + (2 if self.toe_on else 0) \
            + (2 if self.ib_on else 0)
        assert self.n_synapses == sum(expected.values()) == n_syn, \
            f"synapse transcription mismatch: {self.n_synapses} vs {n_syn}"
        assert dict(self.synapse_counts) == expected, "per-tag mismatch"
        assert len(self.idx) == n_neu, f"neuron count {len(self.idx)} != {n_neu}"
        assert len(self.inputs) == n_in, f"port count {len(self.inputs)} != {n_in}"
        assert len(self.muscle_outputs) == 12, "muscle outputs != 12"
        if self.heel_on:
            assert expected.get("heel_rge_in") == 4 and \
                expected.get("heel_pf") == 4
        if self.toe_on:
            assert expected.get("toe_in") == 2 and expected.get("toe_df") == 2 \
                and expected.get("toe_df_mn") == 2
        if self.ib_on:
            assert expected.get("ib_load") == 8
        # the M4 invariants must all still hold underneath
        assert expected.get("pf_drive") == 8 and expected.get("rg_laminate") == 12
        assert expected.get("comm_c1") == 2 and expected.get("c1_SynAmp2.749") == 2


def build(template_path=None, dt: float = DT, gains: dict | None = None,
          tau_rg_nap_h: float | None = None, comm: float = 1.0,
          heel: bool = True, toe: bool = True, ib: bool = True) -> W2LAffNet:
    """Build + compile the afferented split-RG network (Ben's gains by
    default; `gains` overrides individual aff_* / W2L keys)."""
    if template_path is None:
        template_path = SPLIT.SPINAL / "connectome_templates.json"
    tpl = json.loads(Path(template_path).read_text(encoding="utf-8"))
    if "w2laproj" not in tpl:
        raise KeyError(f"{template_path}: no 'w2laproj' key")
    if tau_rg_nap_h is not None:
        SPLIT.NAP["tau_max_h"] = float(tau_rg_nap_h)
    split = SPLIT.make_split_template(tpl["w2laproj"], comm=comm)
    aff = make_aff_template(split, heel=heel, toe=toe, ib=ib)
    net = W2LAffNet(template=aff, comm=comm, gains=gains or {},
                    heel_on=heel, toe_on=toe, ib_on=ib)
    net.compile(dt=dt)
    return net


if __name__ == "__main__":
    import io
    import sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    w = build()
    print(f"neurons: {len(w.idx)}  synapses: {w.n_synapses}  "
          f"inputs: {len(w.inputs)}  outputs: {len(w.muscle_outputs)}")
    print("ports:", w.inputs)
    print("ports inject into:", {p: w.port_targets[p] for p in w.inputs})
    print("synapses per tag:", w.synapse_counts)
    # ablated variants keep the census honest
    w2 = build(heel=False, toe=False, ib=False)
    print(f"all-off: neurons {len(w2.idx)}  synapses {w2.n_synapses}  "
          f"inputs {len(w2.inputs)}  (= M4 coupled 87/166/6)")
