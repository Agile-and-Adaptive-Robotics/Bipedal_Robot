# Modern Robotics nomenclature — dissertation math conventions

Campaign: post-defense revision, 2026-10-15. Target conventions per Lynch & Park,
*Modern Robotics* (the "MR" conventions), as adopted for all dissertation math. Paths
relative to the repo root (`D:\Github\Bipedal_Robot`). Current-state counts below come from
literal searches over `ProofFinal/` run on 2026-10-15: `\Vec` 6 hits (all in
`chapters/20-methods.tex:117-149`), `\mathbf` 46 hits (main.tex macros + 20-methods,
30-results, 92-AppendixA, 94-AppendixC), `\hat` 24 hits (same files), lowercase `\vec`
3 hits (`chapters/92-AppendixA.tex:215,218,243`). The build resolves cleanly — no
"Undefined control sequence" in `ProofFinal/build_wf/main.log` — so `\Vec` is defined by
the loaded math packages (`mathabx` among them, `main.tex:19`); the class defines no `Vec`.

## 1. The rules

1. **Scalars are italic.** $l_0$, $P$, $\theta_{k}$, $r_{k}$, $F$, $M$, $\varepsilon$,
   $\chi_{0}$…$\chi_{3}$.
2. **Vectors are bold lowercase.** $\mathbf{d}$, $\mathbf{p}$, $\mathbf{v}$, and (rulings
   below) $\mathbf{f}$, $\mathbf{m}$. Bold Greek vectors use `\bm{}` or `\boldsymbol{}`
   (`bm` is loaded, `main.tex:18`; `\mathbf` does not bold Greek), e.g. $\bm{\omega}$.
3. **Matrices are bold UPPERCASE.** $\mathbf{R}$, $\mathbf{T}$, $\mathbf{K}$ (bracket
   stiffness), $\mathbf{J}$. Frame transforms use MR subscript pairs:
   $\mathbf{R}_{sb}$, $\mathbf{T}_{sb}$ — target frame replaces both the comma pair
   (`\mathbf{R}_{s,br}`) and the leading-superscript form (`{}^{s}\mathbf{T}_{br}`).
4. **Unit vectors are hatted.** Hat on the bold letter; frame/axis subscripts sit OUTSIDE
   the hat: $\hat{\mathbf{z}}_{k}$, $\hat{\mathbf{x}}_{br}$. Never
   `\hat{z_{k}}` (hat drawn over the subscript — current `20-methods.tex:120`).
5. **Operators are upright.** Already compliant, keep: `\argmin`/`\argmax`
   (`main.tex:61-62`), `\norm{}`/`\abs{}` (`main.tex:121-122`), `\textrm{dist}`
   (`main.tex:120`), built-in trig/log. New operators via `\DeclareMathOperator`, never
   italic. Units stay upright through siunitx (existing practice).
6. **Time derivatives are dots over the symbol.** $\dot{\theta}$, $\dot{l}_\text{m}$,
   $\dot{\mathbf{p}}$ for vectors. Dot notation is already in use
   (`40-discussion.tex:130,135`; `92-AppendixA.tex:26,76`; `93-AppendixB.tex:10,81,105`);
   do not introduce Leibniz $d/dt$ where a dot suffices.
7. **Unchanged (out of MR scope, keep as-is):** the `\mathbb{R}^{n}` family
   (`main.tex:99-108`), calligraphic sets ($\mathcal{N}$, $\mathcal{S}$…), `\mat{}` bmatrix
   shortcut (`main.tex:123`), siunitx quantities.
8. **Only if introduced later** (not currently in the dissertation): MR screws/twists/wrenches
   take calligraphic $\mathcal{S}$, $\mathcal{V}$, $\mathcal{F}$.

## 2. Mapping table — current macros → target forms

| Current use | Where (verified) | Target | Ruling |
|---|---|---|---|
| `\Vec{d}` | 20-methods:117,143,146 | `\mathbf{d}` | distance vector |
| `\Vec{F}` | 20-methods:143,146 | `\mathbf{f}` | vector force goes lowercase per rule 2; the scalar force stays italic $F$ (collision rule below) |
| `\Vec{M}` | 20-methods:143,146 | `\mathbf{m}` | vector moment; scalar torque about $Z$ stays italic $M$ (`20-methods.tex:149-155`) |
| `\Vec{r}` as "moment arm" | 20-methods:117 | italic $r$ (drop vector markup) | the text defines it as "the shortest distance" — a scalar; $r_{k}$ is already scalar at :120-123. This is a semantics fix, not just notation |
| `\Vec{p_{proj,xy}}`, `\Vec{p_f}` | 20-methods:120,123 | `\mathbf{p}_{proj,xy}`, `\mathbf{p}_{f}` | subscript moves outside the bold |
| `\mathbf{R}_z(\theta_z)`, `\mathbf{R}_y(\theta_y)` | 20-methods:181 | keep | already bold uppercase with argument — compliant |
| `\mathbf{R}_{s,br}` | 20-methods:181 | `\mathbf{R}_{sb}` | MR frame subscript pair |
| `{}^{s}\mathbf{T}_{br}` | 20-methods:183 | `\mathbf{T}_{sb}` | leading-superscript frame → subscript pair |
| `\hat{z_{k}}` | 20-methods:120 | `\hat{\mathbf{z}}_{k}` | hat on bold letter, subscript out |
| `\hat z_s`, `\hat x_1`, `\hat y_1`, `\hat y_{br}`, `\hat x_{br}`, `\hat z_{br}` | 20-methods:176-178 | `\hat{\mathbf{z}}_{s}`, `\hat{\mathbf{x}}_{1}`, … | add bold under every unit-vector hat, uniformly |
| `\vec{…}` (3 sites) | 92-AppendixA:215,218,243 | apply the rows above per symbol | eliminate `\vec` from the dissertation |
| `\mathbf{p}_{A}`, `\mathbf{p}_{B}`, `\mathbf{p}_{br}`, `\mathbf{0}^{T}` | 20-methods:178-186, 94-AppendixC passim | keep | already compliant |
| `\mathbf{K}`-style stiffness/compliance blocks | 94-AppendixC:10-155, 92-AppendixA:224-243 | audit only | bold uppercase matrices; check subscripts follow frame-pair style where frames are involved |

Net effect: `\Vec` and `\vec` counts go to zero; every vector is `\mathbf`/`\bm`; every
unit vector is bold with the hat on the letter only.

## 3. Collision rules

- **Bold $\mathbf{f}$ vs. italic $F$:** intentional — bold marks the vector, italic the
  scalar, distinct letters keep the BPA-model equations ($F$, $F_{620}$, $F^{*}$) untouched.
  Within any one equation set, do not also use lowercase italic $f$ for a function; rename
  the function locally ($g$ or $h$).
- **Bold $\mathbf{m}$ vs. mass $m$:** if a derivation needs both, subscript the mass
  ($m_{seg}$-style); one symbol, one meaning per section.
- **Points:** position vectors of named points stay $\mathbf{p}$ with point subscripts
  ($\mathbf{p}_{A}$). MR's bare-italic point style ($q$) is not adopted — the campaign
  rule is minimum churn for attachment-point notation.

## 4. Migration procedure

1. Sweep per block, not globally: `20-methods.tex` §moment-arms/hybrid (:117-155) and the
   bracket-frame block (:176-186 onward), `30-results.tex:257`, `92-AppendixA.tex:182-243`,
   `94-AppendixC.tex:10-155` and :228.
2. Apply the table mechanically. Do not renumber equations; `\label{}`s are unchanged.
3. Recompile through the existing build path; check the log for new errors.
4. Regression greps (literal, over `ProofFinal/chapters/` + `main.tex`):
   - `\Vec` → 0 hits
   - `\vec` → 0 hits
   - `\hat\{[a-zA-Z]+_` → 0 hits (no subscript inside a hat)
   - spot-audit remaining `\mathbf{[A-Z]}` uses for the lowercase-vector rule
5. Report per-block before/after counts in the batch report; Gate 3 of the supervisor
   rubric still applies to any values the touched equations carry.
