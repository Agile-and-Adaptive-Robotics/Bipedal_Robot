# Outline — Appendix A and B drafts

Status: drafted, polished versions in `Notes/overleaf/appendix_A.tex` and
`Notes/overleaf/appendix_B.tex`. Not yet in Overleaf. All content extracted from
the actual MATLAB/AnimatLab sources on 2026-09-07 with file:line provenance.

## Appendix A — Muscle Force and Joint Torque Calculations

Five sections, every equation traceable to code:

1. **BPA force model implementation**
   - Strain measures: eq (app_strain) with chi0, tendon, 2*fitting terms;
     eps620 = (l0-l620)/l0 (code: KMAX); rel strain.
     Provenance: MonoPamDataExplicit.m:210-217, :278, :281; minimizeFlxPin.m:150-164.
   - F* surface (app_fstar) closed form (f_festo.m:1-26); spline/datasheet lookup
     alternative (festo4.m) with clamps; fit provenance bpaFits.m:45-60.
     40 mm coefficients (0.1224, 10.47, 2.023) exist in code only — included,
     flagged as unused in dissertation.
   - Max force arctan (maxBPAforce.m:46/:57); 20 mm CIs (Fmax20.m:26-27);
     fixed 7.5 mm offset; total force product (MonoPamDataExplicit.m:283-287).
2. **Muscle path length and moment arm**
   - Segment-sum path length with frame transforms (:78-95, :122-133).
   - Perpendicular-foot moment arm r = p - u(u·p) (:156-168) and its equivalence
     to the cross-product formula (project_ma) in the text; hypot() scalar
     (predictKneeFlexor20mm.m:61); compressive-region NaN convention (:311).
3. **Compliance model implementation**
   - Corrected length (minimizeFlxPin.m:254, 368-377).
   - **DISCREPANCY — needs Ben's decision:** code uses K_br = diag(chi1, chi2,
     chi2) (:266-268) but dissertation eq (bktstiffness) states diag(chi1, chi2,
     chi1). Polished appendix currently matches the CODE; if Ben decides the
     text is right, flip the polished file's third diagonal entry.
   - Tendon stiffness k = nAE/L with n=2 (x6 for 20 mm), A=1.51e-6 m^2,
     E=193 GPa (:460-471; MonoPam_mult.m:606-617) — not previously in text.
   - fzero force balance on [0, lm-l620], delta=0 guard (:291-305); deflection
     back-update (:324, :334); two-BPA coupled common-force variant
     (MonoPam_mult.m:500-565).
4. **Wrapping-loss implementation**
   - Extensor segment-wise dL with R1=12 mm, R2=40 mm, boundaries 27/-23 deg
     (minimizeExtX3.m Contraction :327-370).
   - 20 mm bend-measure version dL = chi3 * b(theta) * (1-eps*)^2
     (MonoPam_mult.m:398-441; buildKneeFlexorRoute20mm.m:265-268; wRap doc in
     nonlconExclusion.m:150).
   - Loss enters force-producing strain only, not kinematic strain
     (minimizeExtX3.m:236-240).

## Appendix B — Neuron Model Equations

From AnimatLab walker (Biped_2xCPG_wSubs / walk new new tester_Standalone.asim):

1. **Spiking neuron:** leaky integrator membrane eq (app_membrane); spike +
   AHP (1 mV, 3 ms); threshold accommodation (0.3, 10 ms). Parameter table:
   E_rest -60 mV, tau_m 5 ms, V_th -55 mV (motor -59), Ca channel disabled
   (G_max=0).
2. **Non-spiking graded synapses:** clipped-linear conductance (app_nonsyn) +
   current law. Three instances: ACh (E=-10 mV, g=0.5 uS, -65->-20 mV),
   hyperpolarizing IPSP (-70 mV, -55->-40), depolarizing IPSP (-55 mV, -60->-40).
   Electrical synapses rectifying 0.1 / non-rectifying 0.2 available for afferents.
3. **Relation to SNS-Toolbox:** same model classes transfer directly; tonic
   stimulation tests = constant I_app at MN / PF / RG layers (ties to Ch. 6).

## Open decisions for Ben

1. K_br third diagonal entry: fix code, or fix text eq (bktstiffness)? (Polished
   appendix follows the code until you decide.)
2. Include the 40 mm F* coefficients in the appendix? (Currently yes, flagged as
   code-only; harmless, but cuttable.)
3. Appendix B parameter table: confirm these walker values are the ones you want
   documented (they are identical across all neurons per the .asim file).
