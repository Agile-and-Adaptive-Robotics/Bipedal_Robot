import io, re, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
p = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex"
t = open(p, encoding="utf-8").read()
envs = re.findall(r"\\begin\{(\w+\*?)\}", t)
ends = re.findall(r"\\end\{(\w+\*?)\}", t)
from collections import Counter
b, e = Counter(envs), Counter(ends)
print("env mismatch:", {k: (b[k], e[k]) for k in set(b) | set(e) if b[k] != e[k]} or "NONE")
print("equations:", len(re.findall(r"\\begin\{equation\}", t)), "eq-labels:", len(re.findall(r"\\label\{eq:", t)))
dollars = len(re.findall(r"(?<!\\)\$", t))
print("'$' count (must be even):", dollars, "->", "OK" if dollars % 2 == 0 else "ODD!")
keys = ["cut held-out prediction error by two thirds or more",
        "together with the first hardware torque measurements of the redesigned extensor route",
        "where positive", "encoder zero offset, which was removed when the data were loaded",
        "at resting lengths of", "excluded from the pool and the extensor offset allowed negative values",
        "These three tests are the initial-identification set",
        "the shortfall extending slightly past the upper end of the no-tolerance comparison",
        "Because operating at \\qty{620}{\\kPa} across the full desired range of motion",
        "two contact arcs of wrapping radii",
        "scales a physical arc of known size",
        "the selection rule operating as designed rather than an independent confirmation",
        "at or below 0.172",
        "fit this single test more closely still",
        "the zero-correction prediction of the text",
        "an interpretation of the quasistatic balance rather than a quantity this work measures",
        "The surrogate budgets cap run time",
        "whose moving wrap point bends the actuator over a standoff surface",
        "rather than as an isolated pulley effect",
        "when the route builder failed while assembling the redesigned path segments",
        "Investigation, B.B., C.M., and L.B.;",
        "from the direct route's \\qty{44.9}{\\mm} to \\qty{61.2}{\\mm}"]
missing = [k for k in keys if k not in t]
print(f"inserted passages present: {len(keys)-len(missing)}/{len(keys)}; missing: {missing or 'NONE'}")
gone = ["confirmed by first hardware measurements", "halve prediction error on held-out tests across configurations",
        "adjusted so that", "observed plateau", "rigid-body baseline", "on a concatenation error",
        "Investigation, B.B., C.M., and C.M.", "sum of two \\diameter"]
leftover = [k for k in gone if k in t]
print("old wordings fully removed:", "YES" if not leftover else f"LEFTOVER: {leftover}")
print("stray 'at or below 0.17 on' remaining:", "at or below 0.17 on" in t)
