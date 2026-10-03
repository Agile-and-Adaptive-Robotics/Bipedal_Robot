"""Bisect the spiking conversion: neurons-only or synapses-only variant."""
import io
import sys
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

SPIKING_TEMPLATE = (
    "<SynapseType>\n"
    "<Name>{name}</Name>\n"
    "<ID>{sid}</ID>\n"
    "<Type>SpikingChemical</Type>\n"
    "<Equil>{equil}</Equil>\n"
    "<SynAmp>{synamp}</SynAmp>\n"
    "<Decay>10</Decay>\n"
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


def fix_neuron(m):
    block = m.group(0)
    thr_m = re.search(r"<InitialThresh>([-\d.eE+]+)</InitialThresh>", block)
    if not thr_m or float(thr_m.group(1)) < 0:
        return block
    name_m = re.search(r"<Name>([^<]*)</Name>", block)
    name = name_m.group(1) if name_m else ""
    new_thr = -59.0 if re.search(r"\bMN\b", name) else -55.0
    block = re.sub(r"<InitialThresh>[-\d.eE+]+</InitialThresh>",
                   f"<InitialThresh>{new_thr:g}</InitialThresh>", block)
    block = re.sub(r"<RelativeAccom>[^<]*</RelativeAccom>",
                   "<RelativeAccom>0</RelativeAccom>", block)
    return block


def convert(inpath, outpath, do_neurons, do_syn):
    s = open(inpath, encoding="utf-8").read()
    if do_neurons:
        s = re.sub(r"<Neuron>.*?</Neuron>", fix_neuron, s, flags=re.S)
    if do_syn:
        ns_m = re.search(r"<NonSpikingSynapses>(.*?)</NonSpikingSynapses>", s, flags=re.S)
        ns_body = ns_m.group(1)
        converted = []
        for b in re.findall(r"<SynapseType>.*?</SynapseType>", ns_body, flags=re.S):
            name = re.search(r"<Name>([^<]*)</Name>", b).group(1)
            sid = re.search(r"<ID>([^<]*)</ID>", b).group(1)
            equil = re.search(r"<Equil>([^<]*)</Equil>", b).group(1)
            synamp = float(re.search(r"<SynAmp>([^<]*)</SynAmp>", b).group(1))
            converted.append(SPIKING_TEMPLATE.format(
                name=name, sid=sid, equil=equil, synamp=synamp,
                maxrel=max(5.0, 2.0 * synamp)))
        s = s[:ns_m.start(1)] + "\n" + s[ns_m.end(1):]
        sp_close = s.find("</SpikingSynapses>")
        s = s[:sp_close] + "".join(converted) + s[sp_close:]
    with open(outpath, "w", encoding="utf-8", newline="") as f:
        f.write(s)
    print("wrote", outpath, "neurons:", do_neurons, "syn:", do_syn)


if __name__ == "__main__":
    inp, outp, which = sys.argv[1], sys.argv[2], sys.argv[3]
    convert(inp, outp, which in ("n", "both"), which in ("s", "both"))
