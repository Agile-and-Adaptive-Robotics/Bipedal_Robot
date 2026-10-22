import os

p = os.path.join(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation",
    r"Notes\advisor_revision_20261022\patches\README.md",
)
s = open(p, encoding="utf-8").read()
s = s.replace(
    "— standard context diffs. Apply from `ProofFinal/`:",
    "— standard context diffs. Apply from the REPO ROOT (`D:\\Github\\Bipedal_Robot`), quoting the whole path because of the spaces:",
)
s = s.replace(
    "Apply from `ProofFinal/` with (zero-context diffs need the flag):",
    "Apply from the REPO ROOT with (zero-context diffs need the flag):",
)
s = s.replace(
    'git apply --unidiff-zero ".../patches/30-torque/01-physical-validation-sweep.patch"',
    'git apply --unidiff-zero "Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/patches/30-torque/01-physical-validation-sweep.patch"',
)
with open(p, "w", encoding="utf-8", newline="") as f:
    f.write(s)
print("README apply instructions fixed")
