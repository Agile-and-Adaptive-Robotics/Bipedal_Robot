"""Patch Opt_run_Routing27.m: add the too-close soft spacing penalty, reference the
comfort threshold, and adjust the docstring. The hard c5 already uses ctx.spacingRequired
(which is now 1 mm); add a soft penalty for spacing below ctx.spacingComfort."""
from pathlib import Path

p = Path(r"D:\Github\Bipedal_Robot\tmp\routing27_v2_src\Opt_run_Routing27.m")
t = p.read_text(encoding="utf8")

# 1. Add too-close penalty after the JspaceOut block
old_space = """        JspaceOut = ctx.outerW * (max(0, ...
            ctx.spacingRequired - minSpacingOut) / ctx.spacingRequired) .^ 2;
    else
        JspaceOut = 0;
    end"""
new_space = """        JspaceOut = ctx.outerW * (max(0, ...
            ctx.spacingRequired - minSpacingOut) / ctx.spacingRequired) .^ 2;
    else
        JspaceOut = 0;
    end

    % Too-close-to-other-BPAs soft penalty (Ben 2026-10-03: 'penalize if the
    % routes are too close'; the BPAs will squish each other physically, but
    % the optimizer should prefer routes that don't force squishing).
    JspaceClose = 0;
    if nFin > 0 && minSpacing < ctx.spacingComfort
        JspaceClose = 1e-1 * ((ctx.spacingComfort - minSpacing) / ...
            ctx.spacingComfort) .^ 2;
    end"""
assert old_space in t, "spacing block not found"
t = t.replace(old_space, new_space)

# 2. Add JspaceClose to the total objective
old_J = "J = 1e5 * Jworst + 1e3 * Jtorque + Jover + Jreact + Jfoot ..."
new_J = "J = 1e5 * Jworst + 1e3 * Jtorque + Jover + Jreact + Jfoot + JspaceClose ..."
assert old_J in t, "J assembly not found"
t = t.replace(old_J, new_J)

# 3. Add JspaceClose to the result struct
old_res = "'Jrom', stF.Jrom, 'Janat', stF.Janat, ..."
new_res = "'Jrom', stF.Jrom, 'Janat', stF.Janat, 'JspaceClose', stF.JspaceClose, ..."
assert old_res in t, "result struct not found"
t = t.replace(old_res, new_res)

# 4. Add JspaceClose to the st struct
old_st = "'Jrom', Jrom, 'Janat', Janat, ..."
new_st = "'Jrom', Jrom, 'Janat', Janat, 'JspaceClose', JspaceClose, ..."
assert old_st in t, "st struct not found"
t = t.replace(old_st, new_st)

# 5. Update the header docstring
old_doc = "%     - Anatomy / femoral-neck guard: mean over primary-DOF cells of the"
new_doc = ("%     - Anatomy / femoral-neck guard: mean over primary-DOF cells of the\n"
           "%     - Too-close-to-other-BPAs penalty: soft, squared, threshold at the\n"
           "%       BPA-diameter comfort distance (the BPAs will squish each other\n"
           "%       physically; the optimizer should just prefer non-squishing routes).")
assert old_doc in t, "docstring anchor not found"
t = t.replace(old_doc, new_doc)

p.write_text(t, encoding="utf8", newline="\n")
print("Opt_run_Routing27.m patched")
