import io

path = "proposals/30-torque.tex"
with io.open(path, "r", encoding="utf-8-sig", newline="") as f:
    text = f.read()

anchor = "The first physical torque measurements of the redesigned extensor route have been collected."
i = text.find(anchor)
assert i >= 0, "anchor not found"

# find end of that paragraph (next double newline after anchor)
j = text.find("\n\n", i)
assert j >= 0
old_para = text[i:j]

new_para = (
    "Physical validation of the redesigned extensor route is under way with the "
    "$\\diameter$\\qty{20}{\\mm} actuator. In the first test, on the earlier direct "
    "route, three isometric points were held near full extension: measured knee "
    "torque of \\qty{1.163}{\\newton\\cdot\\meter} at a knee angle of \\ang{10}, "
    "\\qty{7.849}{\\newton\\meter} at \\ang{11}, and \\qty{23.437}{\\newton\\meter} "
    "at \\ang{10}, against hung loads of \\qtylist{3.7;25.2;76}{\\newton} at gauge "
    "pressures of \\qtylist{388.5;388;592}{\\kPa}, respectively. The first and third "
    "points sit at essentially the same knee angle and span a twentyfold range of "
    "held torque, so they begin to map the route's response near extension as "
    "pressure and load vary. On the rebuilt cam route, a nine-point isometric sweep "
    "then held loads of \\qtylist{100;116;126;140;168;173.4;180.2;181;181.5}{\\newton} "
    "at knee angles from \\ang{-119} to \\ang{5} and gauge pressures from "
    "\\qty{256}{\\kPa} to \\qty{488}{\\kPa}, each with its load-cell angle recorded "
    "in the tibia frame and its measured wrench carried to the moving knee center by "
    "the same adjoint transform as the identification data, at the measured "
    "reaction-point distance of \\qty{259.9}{\\mm}. The adjoint-reduced torques are "
    "entered for the six most recent holds and span \\qty{28.8}{\\newton\\meter} to "
    "\\qty{50.6}{\\newton\\meter}, with the reduction of the three earliest holds in "
    "progress. The pressures read against the design on their own: the deep-flexion "
    "hold at \\ang{-119} was carried at \\qty{256}{\\kPa}, at the low end of the cam "
    "design's predicted \\qtyrange{269}{507}{\\kPa} requirement band, and the "
    "extension holds at \\ang{-1} and \\ang{5} at \\qty{488}{\\kPa} and "
    "\\qty{480}{\\kPa} against the \\qty{507}{\\kPa} worst case, so the route held "
    "loads across the range without approaching the \\qty{620}{\\kPa} ceiling. The "
    "measured-versus-predicted envelope comparison for the sweep, and the "
    "corresponding flexor measurements, are pending."
)

text = text[:i] + new_para + text[j:]
with io.open(path, "w", encoding="utf-8", newline="") as f:
    f.write(text)
print("replaced paragraph of", len(old_para), "chars with", len(new_para), "chars")
