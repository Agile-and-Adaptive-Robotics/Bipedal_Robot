# Is the additive complement ζ = 1 − ε* justifiable from Chou & Hannaford (1996)?

Working note for the torque paper / dissertation. Question (Ben, 2026-10-18): can the additive
complement term ζ = 1 − ε*, which scales the wrapping loss Δl = χ₃·R·|θ_wrap − θ_k|·ζ²
(eq:deltaL / eq:comp of the manuscript), be related to the length–diameter (L–D) relation of
Chou & Hannaford's braid model — scientifically and mathematically?

## Short answer

**Yes, in a precise and honest sense.** One power of ζ is an exact consequence of the
Chou–Hannaford inextible-thread constraint: the diameter-for-length exchange rate per unit
contraction carries exactly one factor of the remaining-length fraction, dD/dε ∝ (1−ε)/D.
The second power is the standard second-order nature of bending accommodation for a helically
wound, inextensible-thread cylinder (first-order intrados/extrados effects cancel around the
azimuth). Together they motivate Δl ∝ (arc length)·ζ² with an identified coefficient χ₃. The
square is therefore *consistent with* C&H geometry and *parsimonious*, but not uniquely implied
by it — which is exactly why χ₃ is identified from data rather than derived. That is the
defensible claim the paper can carry.

## 1. What Chou & Hannaford give

For a McKibben/BPA braid with n thread turns of inextensible thread length b, muscle length L,
diameter D, and braid angle θ (thread-to-axis):

  b² = L² + (nπD)²,   L = b·cos θ,   nπD = b·sin θ.          (C&H Eq. set)

With contraction ε = (L₀ − L)/L₀ and L = L₀(1 − ε) = b·cos θ₀·(1 − ε):

  cos θ(ε) = (1 − ε)·cos θ₀,   D(ε) = (b/nπ)·√(1 − (1 − ε)²·cos²θ₀).   (exact)

## 2. One power of ζ is exact: the exchange-rate argument

Differentiate the constraint 2L·dL + 2(nπ)²·D·dD = 0 and use dL/dε = −L₀:

  dD/dε = L₀·L / (nπ)²D = L₀²(1 − ε) / (nπ)²D.

So the rate at which the braid can still convert diameter into axial length — the only currency
a contracting BPA has — falls in exact proportion to the remaining-length fraction (1 − ε).
Normalize contraction by the unloaded 620 kPa contraction range ε* = ε/ε_620 as the manuscript
does; the normalized remaining reservoir is ζ = 1 − ε* = (L − L_620)/(L₀ − L_620)
= (cos θ − cos θ_620)/(cos θ₀ − cos θ_620): **ζ is the normalized remaining braid-angle window,
a Chou–Hannaford-native quantity.** Any mechanism that consumes usable contraction (bending
included) must therefore shut down at least linearly in ζ: at full contraction the window is
closed and the loss must vanish. Δl ∝ R·φ·ζ^1 is thus rigorously grounded.

## 3. The second power: second-order bending accommodation

When the actuator wraps over radius R through bend angle φ, points on the extrados would need
arc length (1 + D/2R)·Rφ and points on the intrados (1 − D/2R)·Rφ. A helical thread, however,
spends equal thread length on every azimuth: each turn places half its length on the intrados
and half on the extrados, so the first-order (linear in D/2R) length demands cancel along the
thread, and the leading net effect of curvature on an inextensible helix is second order
(the classical result for wound structures; the same cancellation that makes a helix insensitive
to first-order bending). The braid accommodates the residual second-order demand through the
section-2 exchange, whose capacity itself carries one power of ζ. A demand that is second order
in the accommodation amplitude, met by an exchange whose rate scales as ζ, gives a leading
loss term ∝ ζ². This is an argument of consistency and dimensional parsimony, not a uniqueness
proof: an exponent between 1 and 2 could be constructed from alternative wall-mechanics
assumptions, and the tension-asymmetry mechanism the manuscript describes (intrados wall
unloading under Laplace law on a curved tube) enters the same order.

## 4. What the paper can claim

- ζ = 1 − ε* is the normalized remaining braid-angle window of the C&H braid; through
  dD/dε = L₀²(1 − ε)/(nπ)²D the manuscript's ε_620 normalization is the natural bookkeeping of
  the length–diameter exchange capacity.
- The wrapping loss consequently vanishes at the contraction limit (exact) and its leading
  dependence on the remaining window is second order (second-order azimuthal cancellation of
  bending in an inextible-thread helix).
- χ₃ remains an identified coefficient: C&H fixes the structure (arc length × window²), not the
  magnitude, because the wall thickness, thread compliance, and bladder pressure distribution
  that set the magnitude are outside the ideal C&H model.

Suggested Discussion sentence (journal register, if wanted):
"Expressed in the braid variables of \citet{chou_measurement_1996}, the factor ζ = 1 − ε* is the
fraction of the actuator's contraction window that remains: differentiating the inextible-thread
constraint b² = L² + (nπD)² shows the diameter-for-length exchange rate itself falls in
proportion to the remaining length fraction, so any bending mechanism that consumes usable
contraction must vanish at the contraction limit; and because a helical thread spends equal
length on either side of a bend, the leading accommodation of curvature is second order, which
motivates the quadratic form of Equation~\eqref{eq:deltaL}. The scale of the loss, set by wall
mechanics outside the ideal braid model, is what χ₃ identifies."

## Sources

- Chou, C.-P., & Hannaford, B. (1996). Measurement and modeling of McKibben pneumatic
  artificial muscles. IEEE Transactions on Robotics and Automation, 12(1), 90–102.
  (bib key: chou_measurement_1996 — already in BolenFrontiers22.bib.)
- Manuscript equations: eq:comp (ζ = 1 − ε*), eq:deltaL (Δl = χ₃ R |θ_wrap − θ_k| ζ²),
  eq:x3 (corrected length), ProofFinal/Main.tex.
