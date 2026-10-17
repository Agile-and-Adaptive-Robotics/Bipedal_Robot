import re, io

p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
c = io.open(p, encoding="utf-8").read()

before_phi = len(re.findall(re.escape(r"$\phi$") + r"\s?\\qty", c))
c = re.sub(re.escape(r"$\phi$") + r"\s?\\qty", r"$\\diameter$\\qty", c)
# also bare \phi used as diameter directly attached (e.g. \phi\qty without $)
before_phi2 = len(re.findall(r"\\phi\\qty", c))
c = c.replace(r"\phi\qty", r"\diameter\qty")

# em-dash pair in the two-stage paragraph
old_pair = r"locked $(\chi_{1}, \chi_{2})$ to the flexor front values --- enforcing a single stiffness pair across the flexor and extensor configurations --- and solved"
new_pair = r"locked $(\chi_{1}, \chi_{2})$ to the flexor front values, enforcing a single stiffness pair across the flexor and extensor configurations, and solved"
pair_hit = old_pair in c
c = c.replace(old_pair, new_pair)

# xiProjection caption force-vector notation
old_cap = r"satisfy $\mathbf{F}'_2=-\mathbf{F}'_1$"
new_cap = r"satisfy $\mathbf{f}_2=-\mathbf{f}_1$"
cap_hit = old_cap in c
c = c.replace(old_cap, new_cap)

io.open(p, "w", encoding="utf-8").write(c)
print("phi$..qty replaced:", before_phi, "| bare phi-qty replaced:", before_phi2)
print("em-dash pair fixed:", pair_hit, "| xiProjection caption fixed:", cap_hit)
print("remaining $\\phi$:", c.count(r"$\phi$"), "| remaining ---:", c.count("---"),
      "| remaining F-prime:", c.count(r"\mathbf{F}'"))
