"""Convert a non-spiking AnimatLab .asim to a spiking copy (_spiking suffix).

Mechanism (verified against the AnimatLab public source: IntegrateFireSim has ONE
serialized Neuron class; spiking is dynamic threshold crossing, m_bSpike=(Vm>Thresh)):
  - every <Neuron> whose InitialThresh >= 0 (unreachable: 50/200) becomes spiking:
      InitialThresh -> -55 (Li-native) or -59 for MN-named cells (Li MN threshold)
      RelativeAccom -> 0   (Li's spiking cells turn accommodation off)
    everything else (rest, tc, GMaxCa, AHP, tonic, Ca gates, IDs) unchanged.
  - every NonSpikingChemical <SynapseType> becomes SpikingChemical:
      keeps Name/ID/Equil/SynAmp (SynAmp = per-spike conductance increment, uS)
      adds Decay=10 ms, MaxRelCond=max(5, 2*SynAmp), neutral Hebbian/facilitation,
      VoltDep off; block moves into <SpikingSynapses>; <NonSpikingSynapses> left empty.
  - connexions keep their SynapseTypeID references (IDs unchanged) so the wiring
    and all chart-column GUIDs stay valid.
Electrical synapse types and ExternalStimuli untouched.

Usage: python make_spiking.py <in.asim> <out.asim>
"""
import io
import sys
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

_converted_ids = []

# v2 calibrated thresholds: baseline graded peak - 1.5 mV per cell class,
# peaks measured from the baseline runs (baseline_graded_ranges.json):
#   RG/PF/IN/Ia/Ib peaks ~ -52 mV  ->  -53.5
#   MN peaks ~ -42.4 mV            ->  -44
#   RE (Renshaw) peaks ~ -43 mV    ->  -44.5
THR_V1 = {"MN": -55.0, "default": -55.0}  # v1 actually used -59 for MN, -55 else
CAL_THR = {"MN": -44.0, "RE": -44.5, "default": -53.5}

SPIKING_TEMPLATE = (
    "<SynapseType>\n"
    "<Name>{name}</Name>\n"
    "<ID>{sid}</ID>\n"
    "<Type>SpikingChemical</Type>\n"
    "<Equil>{equil}</Equil>\n"
    "<SynAmp>{synamp}</SynAmp>\n"
    "<Decay>{decay}</Decay>\n"
    "<RelFacil>1</RelFacil>\n"
    "<FacilDecay>100</FacilDecay>\n"
    "<VoltDep>False</VoltDep>\n"
    "<MaxRelCond>{maxrel}</MaxRelCond>\n"
    "<SatPSPot>-30</SatPSPot>\n"
    "<ThreshPSPot>-60</ThreshPSPot>\n"
    "<Hebbian>False</Hebbian>\n"
    "<MaxAugCond>1</MaxAugCond>\n"
    "<LearningInc>0.1</LearningInc>\n"
    "<LearningTime>20</LearningTime>\n"
    "<AllowForget>False</AllowForget>\n"
    "<ForgetTime>10000</ForgetTime>\n"
    "<Consolidation>1</Consolidation>\n"
    "</SynapseType>\n"
)


def convert(inpath, outpath, calibrated=False, decay=10):
    s = open(inpath, encoding="utf-8").read()
    stats = {"neurons_converted": 0, "mn_neurons": 0, "types_converted": 0,
             "synamp_values": [], "equil_values": [], "calibrated": calibrated}
    _converted_ids.append(set())

    # ---------- neurons ----------
    def fix_neuron(m):
        block = m.group(0)
        name_m = re.search(r"<Name>([^<]*)</Name>", block)
        thr_m = re.search(r"<InitialThresh>([-\d.eE+]+)</InitialThresh>", block)
        if not thr_m:
            return block
        thr = float(thr_m.group(1))
        if thr < 0:  # already spiking (Li MV cells etc.)
            return block
        name = name_m.group(1) if name_m else ""
        if calibrated:
            if re.search(r"\bMN\b", name):
                new_thr = CAL_THR["MN"]
            elif re.search(r"\bRE\b", name):
                new_thr = CAL_THR["RE"]
            else:
                new_thr = CAL_THR["default"]
        else:
            new_thr = -59.0 if re.search(r"\bMN\b", name) else -55.0
        if re.search(r"\bMN\b", name):
            stats["mn_neurons"] += 1
        block = re.sub(r"<InitialThresh>[-\d.eE+]+</InitialThresh>",
                       f"<InitialThresh>{new_thr:g}</InitialThresh>", block)
        block = re.sub(r"<RelativeAccom>[^<]*</RelativeAccom>",
                       "<RelativeAccom>0</RelativeAccom>", block)
        stats["neurons_converted"] += 1
        return block

    s = re.sub(r"<Neuron>.*?</Neuron>", fix_neuron, s, flags=re.S)

    # ---------- synapse types ----------
    ns_m = re.search(r"<NonSpikingSynapses>(.*?)</NonSpikingSynapses>", s, flags=re.S)
    if not ns_m:
        raise SystemExit("no NonSpikingSynapses section found")
    ns_body = ns_m.group(1)
    blocks = re.findall(r"<SynapseType>.*?</SynapseType>", ns_body, flags=re.S)
    converted = []
    for b in blocks:
        name = re.search(r"<Name>([^<]*)</Name>", b).group(1)
        sid = re.search(r"<ID>([^<]*)</ID>", b).group(1)
        equil = re.search(r"<Equil>([^<]*)</Equil>", b).group(1)
        synamp = float(re.search(r"<SynAmp>([^<]*)</SynAmp>", b).group(1))
        maxrel = max(5.0, 2.0 * synamp)
        _converted_ids[0].add(sid)
        converted.append(SPIKING_TEMPLATE.format(
            name=name, sid=sid, equil=equil, synamp=synamp, maxrel=maxrel, decay=decay))
        stats["types_converted"] += 1
        stats["synamp_values"].append(synamp)
        stats["equil_values"].append(equil)

    # empty the NonSpikingSynapses section
    s = s[:ns_m.start(1)] + "\n" + s[ns_m.end(1):]
    # splice converted blocks into SpikingSynapses before its close
    sp_close = s.find("</SpikingSynapses>")
    if sp_close < 0:
        raise SystemExit("no SpikingSynapses section found")
    s = s[:sp_close] + "".join(converted) + s[sp_close:]

    # ---------- connexions: flip their own kind enum <Type>1</Type> (NonSpikingChemical)
    # to <Type>0</Type> (SpikingChemical) for every connexion referencing a converted type
    converted_ids = {b_id for b_id in _converted_ids[0]}
    flipped = [0]

    def fix_connexion(m):
        block = m.group(0)
        sid_m = re.search(r"<SynapseTypeID>([^<]*)</SynapseTypeID>", block)
        if sid_m and sid_m.group(1) in converted_ids:
            new = re.sub(r"<Type>1</Type>", "<Type>0</Type>", block)
            if new != block:
                flipped[0] += 1
            return new
        return block

    s = re.sub(r"<Connexion>.*?</Connexion>", fix_connexion, s, flags=re.S)
    stats["connexions_flipped"] = flipped[0]

    with open(outpath, "w", encoding="utf-8", newline="") as f:
        f.write(s)
    return stats


if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    cal = "--calibrated" in sys.argv
    dec = 10
    if "--decay" in sys.argv:
        dec = float(sys.argv[sys.argv.index("--decay") + 1])
    st = convert(inp, outp, calibrated=cal, decay=dec)
    print(f"wrote {outp}" + (" [calibrated v2]" if cal else " [v1]"))
    print(f"  neurons converted: {st['neurons_converted']} (MN@-59: {st['mn_neurons']})")
    print(f"  synapse types converted: {st['types_converted']}")
    print(f"  connexions flipped to Type=0: {st.get('connexions_flipped', 0)}")
    print(f"  SynAmp range kept: {min(st['synamp_values'])} .. {max(st['synamp_values'])}")
