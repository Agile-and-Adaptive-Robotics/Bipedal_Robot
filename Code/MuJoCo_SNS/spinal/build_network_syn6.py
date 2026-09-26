"""VARIANT 2 of the s3k walker architecture: the 6-synergy network (syn6).

Campaign 2026-09-25 (goal 4). Selected by env ``AARL_NET=syn6`` OR
``params.G["syn6"] > 0`` (default 0 = the stock build_network.build runs;
existing topology untouched at defaults - regression-gated at
410 neurons / 376 inputs / 1186 synapses).

Architecture (per side unless noted):

  ONE rhythm generator, identical to the stock ``_build_rg`` pattern:
  persistent-Na half-centers RG_E/RG_F (params.NAP, fixed tau_h), IN-
  laminated mutual inhibition through InE/InF, DRIVE/POSTURE edges.

  SIX pattern-formation layers PF_S1..PF_S6, one per NMF synergy channel
  of ``synergy_basis.npz`` (W_r/W_l [43 muscles x 6], name-keyed). Each
  E-family channel is driven by RG_E, each F-family channel by RG_F.
  Family membership and tau staggering come from the measured per-synergy
  phase profiles (fsa_results/fsa_backsolve.npz ``{side}_pf_phase_mean``,
  heel-strike=0 / toe-off=50 convention): stance-dominant channels are
  E-family, swing-dominant F-family; within a family the channels are
  ranked by peak phase and get tau multipliers interpolating the stock
  PF_SHAPE staggering (E: 0.6/0.9/1.2, F: 0.9/1.2/1.6). Reciprocal PF
  lamination: E channels -> PF_IN_E -> (inh) F channels and back, same
  topology as the stock four-cell PF.

  PF -> MN conductances are DATA, not searched: for synergy k and muscle
  act, g = Eq18(W[act, k]) with the Szczecinski 2017 Eq-18 mapping
  implemented and audited in fsa_backsolve.analytical_conductance
  (g = k*R*Gm/(dE - k*R), R = E_HI = 5 mV, Gm = 1 uS, dE = E_REV_EXC =
  8 mV; imported lazily from fsa_backsolve so the math cannot drift).
  Edges exist only where W > 0 (conditional topology). All 86 non-pruned
  W entries have k < dE/R = 1.6, so no entry is Eq-18-invalid. The 6
  runner-pruned actuators (runner.PRUNE_MUSCLES: quad_fem/gem/peri r+l)
  keep their MN + afferent arc (the runner indexes mn_names for every
  actuator) but receive NO synergy drive: their names are absent from the
  W basis, exactly mirroring the runner's force-zeroing of the same list.

  Rule-file dress (gains verbatim from Ben's connectome JSONs unless a
  deviation is noted):
  - HEEL = stance-phase reset AT the PF layer (master rules
    ben_rules_20260924.json): HEEL_IN -> InE exc 0.5, -> InF exc 0.5,
    -> PF_IN_E exc 0.5. TOE = dorsiflexion inhibition ONLY:
    TOE_IN -> TOEDF exc 5, TOEDF -> the dorsiflexion channel (the
    F-family synergy with the largest ankle_df W mass) inh 2.749.
  - Ib load (master rules): per stance-group muscle Ib -> its group
    IBEXC IN exc 0.5, IBEXC -> homonymous MN exc 0.5 (stance reversal);
    per side the LOAD IN (port LOAD_c) collects stance Ib (Ib -> LBIN
    exc 0.5) and drives RG_E exc 0.5, InE exc 0.5 and each E-family PF
    channel exc 0.5. (Drawing: "Ib grp -> RG-E 0.5, -> IN-InE 0.5,
    -> HC-PF-E 0.5".)
  - Shevtsova commissurals (ben_shevtsova_20260924.json, symmetric):
    RG-F -> V2a exc 1.0, RG-F -> V0D exc 0.7, RG-E -> V3-E exc 0.35,
    V2a -> V0V exc 1.0, V0V -> contra InE exc 0.6 (Ini folded into the
    existing laminated InE, deviation D2), V0D -> contra RG-F inh 0.07,
    V3-E -> contra RG-E exc 0.02 (the file is asymmetric 0.02/0.5; the
    task text's "small exc" reading is used, deviation D2), V3-E ->
    contra InE exc 1.0. Brainstem gamma/alpha are built as shared
    neurons driven by DRIVE (the runner's only descending port -
    deviation D3): DRIVE -> GAMMA/ALPHA exc 1.0, GAMMA -> RG-E exc 0.5,
    ALPHA -> RG-F exc 0.5, ALPHA -> V0V/V0D inh 0.5.
  - Shinohara autogenic afferent expansion (master rules motif, per
    muscle): Ia -> MN homo exc 2.0; Ia -> IaIN exc 1.0, IaIN -> antagonist
    MNs inh 0.5, IaIN <-> antagonist IaIN inh 0.5; II -> IIX exc 1.0,
    IIX -> MN homo exc 0.5; II -> IIIN exc 0.5, IIIN -> antagonist MNs
    inh 0.5; Ib -> IBIN exc 1.0, IBIN -> MN homo inh 0.5, IBIN <->
    antagonist IBIN inh 0.5. Proprioception feeds back to ITSELF
    (autogenic, per muscle) - never group-broadcast (Ben's reading rule).

  Deviations from the rule-file numbers (all documented, none silent):
  D1: RG->PF and PF reciprocal lamination use the stock S3K-tuned knobs
      G["rg_to_pf"] (2.4) / G["pf_recip_inh"] (4.0) instead of the
      drawing's 0.1 / 2.749 - those are W2L spiking-template scale; the
      Eq-18 mapping presumes the PF source itself reaches E_HI, which is
      the S3K drive scale (fsa analytic_pred peak 1.27 confirms).
  D2: Shevtsova IniE/IniF/Ini folded into the existing laminated InE/InF
      (the RG-F->IniF->RG-E and RG-E->IniE->RG-F edges duplicate the
      lamination at the stock G["rg_mutual_inh"] strength).
  D3: brainstem gamma/alpha are DRIVE-driven shared cells (the runner
      exposes one descending port, "DRIVE").

  Runner compatibility: the input-port list matches the stock builder
  for every port the runner writes (DRIVE, POSTURE, BAL_PF, BAL_DF,
  BAL_TRK_EXT/FLX, BAL_LAT_R/L, per side HEEL_c/TOE_c/LOAD_c, per
  actuator POST/Ia/II/Ib). stance_fb=True (mechanosensor INs exist),
  aff_loops/vest False (those pathways are not in this variant's dress).
  net.idx exposes RG_E_{s}, RG_F_{s}, MN_{act} and PF_S{k}_{s}.

  Not implemented (out of this variant's scope): Renshaw cells (not in
  the goal-4 dress list; the stock G["renshaw"] key stays default-0 and
  stock-build-only), AFF_E/AFF_F semi-closed loops, VEST cells, F1/KINH
  swing suppression (no phase-window cells here).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import build_network as SNSpinal
from build_network import TAU, _syn
from muscle_map import EXTENSOR_STANCE_GROUPS, MuscleInfo, classify
from params import G, NAP

from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.networks import Network

HERE = Path(__file__).parent

# Six synergy channels (basis of record; synergy_model.py K = 6).
N_SYN = 6

# Verbatim runner.PRUNE_MUSCLES (runner.py:70-71); not imported to avoid
# pulling the runner/mujoco chain into the builder. Equality with the
# runner's set is asserted by the gate script, not assumed.
PRUNE_MUSCLES = {"quad_fem_r", "quad_fem_l", "gem_r", "gem_l",
                 "peri_r", "peri_l"}

SYNERGY_NPZ = HERE / "synergy_basis.npz"
FSA_PHASE_NPZ = HERE / "fsa_results" / "fsa_backsolve.npz"

# tau multipliers within a family, ranked by measured peak phase
# (early -> late). Interpolated if a family holds fewer/more channels.
SYN6_TAU_E = (0.6, 0.9, 1.2)
SYN6_TAU_F = (0.9, 1.2, 1.6)

# stance-fraction band inside which a channel counts as genuinely MIXED
# and receives BOTH RG drives (soft membership). Outside the band the
# channel is single-family: a small opposite-family admixture on the
# pure channels keeps them faintly on through the antagonistic PF
# lamination, which pins the other family below threshold (the same
# co-contraction deadlock the stock net hit with its old DRIVE->PF
# tonic term; measured 2026-09-25 gate c1 first attempt).
SYN6_MIXED_BAND = (0.35, 0.65)

# Rule-file dress gains (verbatim; see module docstring for D1-D3).
DRESS = dict(
    heel_to_ine=0.5, heel_to_inf=0.5, heel_to_pf_in_e=0.5,
    toe_to_indf=5.0, df_inh=2.749,
    ib_to_ibexc=0.5, ibexc_to_mn=0.5,
    ib_to_lbin=0.5, lbin_to_rge=0.5, lbin_to_ine=0.5, lbin_to_pf_e=0.5,
    ia_to_mn=2.0, ia_to_iain=1.0, iain_to_ant_mn=0.5, iain_mut=0.5,
    ii_to_iix=1.0, iix_to_mn=0.5,
    ii_to_iiin=0.5, iiin_to_ant_mn=0.5,
    ib_to_ibin=1.0, ibin_to_mn=0.5, ibin_mut=0.5,
    rgf_to_v2a=1.0, rgf_to_v0d=0.7, rge_to_v3=0.35,
    v2a_to_v0v=1.0, v0v_to_contra_ine=0.6,
    v0d_to_contra_rgf=0.07,
    v3_to_contra_rge=0.02, v3_to_contra_ine=1.0,
    drive_to_bs=1.0, gamma_to_rge=0.5, alpha_to_rgf=0.5,
    alpha_to_v=0.5,
)

_cached_tables: dict | None = None


def _load_synergy_tables() -> dict:
    """W matrices name-keyed per side + family/tau/df-channel assignment.

    Family: stance-dominant -> 'E', swing-dominant -> 'F', decided from
    the measured phase profiles (fsa phase mean, heel=0/toe-off=50),
    SYMMETRIZED across sides (the mean of the two sides' stance
    fractions). This is the documented mitigation for the known
    S5/S6 bilateral phase instability (synergy_model.py NOTE open item
    (b), measured here: stance_frac S5 = 0.58 r vs 0.24 l): both legs
    run the SAME family assignment, so the rhythm sees one consistent
    pattern; S5 lands F (mean 0.41). Fallback if the fsa file is
    unavailable: per-side W stance/swing mass (asymmetric, flagged).
    tau multiplier: within-family rank by that side's measured peak
    phase, interpolated over SYN6_TAU_E / SYN6_TAU_F. df-channel: the
    side's F channel with the largest ankle_df W mass (the
    dorsiflexion HC the TOE edge inhibits).
    """
    global _cached_tables
    if _cached_tables is not None:
        return _cached_tables
    data = np.load(SYNERGY_NPZ, allow_pickle=True)
    tables: dict[str, dict] = {}
    phase: dict[str, dict] = {}
    for side in ("r", "l"):
        names = [str(x) for x in data[f"muscle_names_{side}"]]
        W = np.asarray(data[f"W_{side}"], dtype=float)
        assert W.shape == (len(names), N_SYN), \
            f"W_{side} shape {W.shape} != ({len(names)}, {N_SYN})"
        # per-synergy stance/swing/df mass over the 43 named muscles
        import muscle_map as mm
        stance_mass = np.zeros(N_SYN)
        swing_mass = np.zeros(N_SYN)
        df_mass = np.zeros(N_SYN)
        group_mass: dict[str, float] = {}
        for i, act in enumerate(names):
            info = classify(act)
            if info is None:
                raise ValueError(f"synergy basis name {act!r} not classifiable")
            prim = info.groups[0]
            group_mass[prim] = group_mass.get(prim, 0.0) + W[i].sum()
            for k in range(N_SYN):
                w = float(W[i, k])
                if w <= 0.0:
                    continue
                if prim in mm.EXTENSOR_STANCE_GROUPS:
                    stance_mass[k] += w
                else:
                    swing_mass[k] += w
                if prim == "ankle_df":
                    df_mass[k] += w
        tables[side] = dict(names=names, W=W, stance_mass=stance_mass,
                            swing_mass=swing_mass, df_mass=df_mass,
                            group_mass=group_mass)
        # measured phase profiles when available (heel=0 / toe-off=50)
        try:
            fsa = np.load(FSA_PHASE_NPZ, allow_pickle=True)
            pm = np.asarray(fsa[f"{side}_pf_phase_mean"], dtype=float)
            # NOTE the grid key has no "pf" in the fsa npz (written as
            # f"{side}_phase_grid", fsa_backsolve.py save block) - using
            # the wrong key here would silently fall back to W masses.
            grid = np.asarray(fsa[f"{side}_phase_grid"], dtype=float)
            assert pm.shape[0] == grid.size and pm.shape[1] == N_SYN
            phase[side] = dict(
                stance_frac=pm[grid < 50].mean(axis=0) / np.maximum(
                    pm[grid < 50].mean(axis=0)
                    + pm[grid >= 50].mean(axis=0), 1e-12),
                peak_phase=np.array([float(grid[int(np.argmax(pm[:, k]))])
                                     for k in range(N_SYN)]))
        except (FileNotFoundError, KeyError, ValueError, AssertionError):
            phase[side] = {}
    # ---- families: symmetrized across sides when both measured ----
    if all(phase.get(s) for s in ("r", "l")):
        mean_sf = 0.5 * (phase["r"]["stance_frac"]
                         + phase["l"]["stance_frac"])
    else:
        mean_sf = None
    for side in ("r", "l"):
        tab = tables[side]
        if mean_sf is not None:
            sf = mean_sf
        else:  # fallback: this side's W stance mass (flagged asymmetric)
            sf = tab["stance_mass"] / np.maximum(
                tab["stance_mass"] + tab["swing_mass"], 1e-12)
        fam = ["E" if sf[k] > 0.5 else "F" for k in range(N_SYN)]
        ph = phase.get(side) or {}
        peak_phase = ph.get("peak_phase",
                            np.arange(N_SYN, dtype=float) * 10.0)
        tau_mult = np.zeros(N_SYN)
        for family, taus in (("E", SYN6_TAU_E), ("F", SYN6_TAU_F)):
            idxs = [k for k in range(N_SYN) if fam[k] == family]
            if not idxs:
                continue
            idxs.sort(key=lambda k: peak_phase[k])
            m = len(idxs)
            for pos, k in enumerate(idxs):
                frac = pos / max(m - 1, 1)
                tau_mult[k] = float(taus[0] + (taus[-1] - taus[0]) * frac)
        f_chans = [k for k in range(N_SYN) if fam[k] == "F"]
        df_channel = (int(f_chans[int(np.argmax(tab["df_mass"][f_chans]))])
                      if f_chans and tab["df_mass"][f_chans].sum() > 0
                      else None)
        tab.update(family=fam, tau=tau_mult, stance_frac=sf,
                   peak_phase=peak_phase, df_channel=df_channel,
                   symmetrized=mean_sf is not None)
    _cached_tables = tables
    return tables


def _eq18_conductances(k_vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Szczecinski 2017 Eq-18, imported from the audited implementation
    (fsa_backsolve.analytical_conductance) so the math cannot drift."""
    from fsa_backsolve import analytical_conductance  # lazy (sklearn heavy)
    return analytical_conductance(k_vals)


@dataclass
class Syn6Network(SNSpinal.SpinalNetwork):
    """Drop-in replacement for SpinalNetwork under AARL_NET=syn6."""

    # evidence attributes (read by the gate scripts / report)
    syn_families: dict = field(default_factory=dict)
    syn_tau_mult: dict = field(default_factory=dict)
    syn_df_channel: dict = field(default_factory=dict)
    syn_g_stats: dict = field(default_factory=dict)
    syn_mixed: dict = field(default_factory=dict)
    pf_cells: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ build
    def __post_init__(self):
        self.net = Network(name="gait2392 spinal syn6")
        n = self.net
        self.f1_kneext_inh = False
        self.phase_reset = False
        self.renshaw = bool(G["renshaw"] > 0.0)   # default 0: absent
        self.stance_fb = True      # HEEL/TOE/LOAD ports + dress (always on)
        self.ia_in = False         # motif is built unconditionally instead
        self.aff_loops = False
        self.vest = False
        self.joint_pf = False

        tables = _load_synergy_tables()

        # ---- descending / balance cells (shared; stock port names) ----
        for name in ("DRIVE", "POSTURE", "BAL_PF", "BAL_DF",
                     "BAL_TRK_EXT", "BAL_TRK_FLX", "BAL_LAT_R", "BAL_LAT_L"):
            self._add(name, TAU["descend"], n)
            n.add_input(name)
            self.inputs.append(name)
        # brainstem gamma/alpha (Shevtsova file; DRIVE-driven, D3) —
        # OFF by default (G["syn6_brainstem"] = 0, JSON-rule key): the
        # runner's single DRIVE port already drives RG-E/RG-F via
        # descend_to_rg_e/f at the S3K-tuned 4.0 nA walk level, and the
        # folded +0.5 gamma/alpha bonus on top measurably E-latched the
        # network in the full runner air run (gate c2 first attempt:
        # E-PF saturated 5.03 mV tonic, RG swing 0.43 mV). Built only
        # when the gain is set (conditional topology).
        self.brainstem = bool(G.get("syn6_brainstem", 0.0) > 0.0)
        if self.brainstem:
            for name in ("BS_GAMMA", "BS_ALPHA"):
                self._add(name, TAU["descend"], n)
            n.add_connection(_syn(DRESS["drive_to_bs"], exc=True),
                             "DRIVE", "BS_GAMMA")
            n.add_connection(_syn(DRESS["drive_to_bs"], exc=True),
                             "DRIVE", "BS_ALPHA")

        # ---- per-side circuitry ----
        for side in self.sides:
            self._build_rg(n, side)
            self._build_pf_s6(n, side, tables[side])
            self._build_commissurals_side(n, side)
        # ---- crossed commissural wiring AFTER both sides exist (the
        # stock builder's ordering rule; sides are built 'l' first so
        # contra targets do not exist yet inside the per-side loop) ----
        if self.interleg:
            for a, b in (("l", "r"), ("r", "l")):
                n.add_connection(
                    _syn(DRESS["v0v_to_contra_ine"], exc=True),
                    f"V0V_{a}", f"InE_{b}")
                n.add_connection(
                    _syn(DRESS["v0d_to_contra_rgf"], exc=False),
                    f"V0D_{a}", f"RG_F_{b}")
                n.add_connection(
                    _syn(DRESS["v3_to_contra_rge"], exc=True),
                    f"V3E_{a}", f"RG_E_{b}")
                n.add_connection(
                    _syn(DRESS["v3_to_contra_ine"], exc=True),
                    f"V3E_{a}", f"InE_{b}")
            # brainstem alpha gates the crossed V0V/V0D pathway (inh 0.5)
            # (built only when G["syn6_brainstem"] > 0)
            if self.brainstem:
                for side in self.sides:
                    n.add_connection(_syn(DRESS["alpha_to_v"], exc=False),
                                     "BS_ALPHA", f"V0V_{side}")
                    n.add_connection(_syn(DRESS["alpha_to_v"], exc=False),
                                     "BS_ALPHA", f"V0D_{side}")

        # ---- muscles: neurons first, then wire ----
        self._order = {a: i for i, a in enumerate(self.muscles)}
        for act, mi in self.muscles.items():
            self._add_muscle_neurons(n, act, mi)
        for act, mi in self.muscles.items():
            self._wire_muscle_s6(n, act, mi, tables.get(mi.side))
        self._wire_balance(n)

    # ------------------------------------------------------------------ parts
    def _build_rg(self, n: Network, side: str):
        """Same pattern as stock _build_rg (persistent-Na half-centers,
        laminated mutual inhibition, DRIVE/POSTURE edges). The v11
        mechanosensor INs are built in the dress section below with the
        master-rules gains instead of the stock gain-keyed edges."""
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        ine, inf = f"InE_{side}", f"InF_{side}"
        nap = NonSpikingNeuronWithPersistentSodiumChannel(
            membrane_capacitance=TAU["rg"], membrane_conductance=1.0,
            resting_potential=0.0, bias=0.0,
            g_ion=np.array([NAP["g_ion"]]),
            e_ion=np.array([NAP["e_ion"]]),
            k_m=np.array([NAP["k_m"]]),
            slope_m=np.array([NAP["slope_m"]]),
            e_m=np.array([NAP["e_m"]]),
            k_h=np.array([NAP["k_h"]]),
            slope_h=np.array([NAP["slope_h"]]),
            e_h=np.array([NAP["e_h"]]),
            tau_max_h=np.array([TAU["rg_nap_h"]]))
        n.add_neuron(nap, name=rg_e)
        self.idx[rg_e] = len(self.idx)
        n.add_neuron(nap, name=rg_f)
        self.idx[rg_f] = len(self.idx)
        for name in (ine, inf):
            self._add(name, TAU["rg"], n)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_e, ine)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), ine, rg_f)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=True), rg_f, inf)
        n.add_connection(_syn(G["rg_mutual_inh"], exc=False), inf, rg_e)
        if G["rg_weak_exc"] > 0.0:
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_e, rg_f)
            n.add_connection(_syn(G["rg_weak_exc"], exc=True), rg_f, rg_e)
        n.add_connection(_syn(G["descend_to_rg_e"], exc=True), "DRIVE", rg_e)
        n.add_connection(_syn(G["descend_to_rg_f"], exc=True), "DRIVE", rg_f)
        n.add_connection(_syn(G["posture_to_rg_e"], exc=True),
                         "POSTURE", rg_e)
        # brainstem dress (Shevtsova): gamma -> RG-E, alpha -> RG-F
        # (built only when G["syn6_brainstem"] > 0; see __post_init__)
        if self.brainstem:
            n.add_connection(_syn(DRESS["gamma_to_rge"], exc=True),
                             "BS_GAMMA", rg_e)
            n.add_connection(_syn(DRESS["alpha_to_rgf"], exc=True),
                             "BS_ALPHA", rg_f)

        # ---- mechanosensor INs (ports HEEL_c/TOE_c/LOAD_c; the runner
        # writes them whenever stance_fb, scaled by G heel_rge/toe_rge/
        # ib_rge which stay default-0 tunables per the JSON RULE) ----
        heel_in, toe_in, lbin = f"HEEL_{side}", f"TOE_{side}", f"LBIN_{side}"
        self._add(heel_in, TAU["preset"], n)
        self._add(toe_in, TAU["preset"], n)
        self._add(lbin, TAU["ib_exc"], n)
        n.add_input(heel_in)
        self.inputs.append("HEEL_c_" + side)
        n.add_input(toe_in)
        self.inputs.append("TOE_c_" + side)
        n.add_input(lbin)
        self.inputs.append("LOAD_c_" + side)

    def _build_pf_s6(self, n: Network, side: str, tab: dict):
        """Six synergy PF layers + reciprocal lamination + the heel/toe
        dress at the PF layer (master rules)."""
        fam, tau_mult = tab["family"], tab["tau"]
        cells = []
        for k in range(N_SYN):
            pf = f"PF_S{k + 1}_{side}"
            self._add(pf, TAU["pf"] * float(tau_mult[k]), n)
            cells.append(pf)
        self.pf_cells[side] = cells
        self.syn_families[side] = list(fam)
        self.syn_tau_mult[side] = [float(t) for t in tau_mult]
        self.syn_df_channel[side] = tab["df_channel"]
        rg_src = {"E": f"RG_E_{side}", "F": f"RG_F_{side}"}
        # RG drive split by each channel's MEASURED stance fraction, but
        # only for genuinely MIXED channels (SYN6_MIXED_BAND): pure
        # channels get the single-family edge. A blind split put a faint
        # opposite-family tone on every channel, which held the
        # antagonistic PF lamination on and pinned the F family below
        # threshold (gate c1 first attempt). The mixed channel - the
        # S5/S6-unstable one - gets a genuine admixture so its waveform
        # is distinct from its family prototypes (pre-fix within-family
        # corr +0.98..+0.997).
        sf = tab["stance_frac"]
        for k in range(N_SYN):
            s = float(sf[k])
            if SYN6_MIXED_BAND[0] < s < SYN6_MIXED_BAND[1]:
                n.add_connection(_syn(G["rg_to_pf"] * s, exc=True),
                                 rg_src["E"], cells[k])
                n.add_connection(_syn(G["rg_to_pf"] * (1.0 - s), exc=True),
                                 rg_src["F"], cells[k])
            elif s >= 0.5:
                n.add_connection(_syn(G["rg_to_pf"], exc=True),
                                 rg_src["E"], cells[k])
            else:
                n.add_connection(_syn(G["rg_to_pf"], exc=True),
                                 rg_src["F"], cells[k])
        # PF lamination (stock topology, generalized to the COMMITTED
        # cells only). A genuinely MIXED channel (SYN6_MIXED_BAND) is
        # excluded: riding both RG drives it would be on through BOTH
        # phases, hold its family IN latched and pin the other family
        # below threshold (gate c1 second attempt measured exactly
        # that). Excluded, it runs as a distinct mixed waveform that
        # cannot latch the layer.
        mixed = [(SYN6_MIXED_BAND[0] < float(sf[k]) < SYN6_MIXED_BAND[1])
                 for k in range(N_SYN)]
        self.syn_mixed[side] = mixed
        pf_in_e = f"PF_IN_E_{side}"
        pf_in_f = f"PF_IN_F_{side}"
        self._add(pf_in_e, TAU["pf"], n)
        self._add(pf_in_f, TAU["pf"], n)
        for k in range(N_SYN):
            if mixed[k]:
                continue
            n.add_connection(_syn(G["pf_recip_inh"], exc=True), cells[k],
                             pf_in_e if fam[k] == "E" else pf_in_f)
        for k in range(N_SYN):
            if mixed[k]:
                continue
            src = pf_in_e if fam[k] == "F" else pf_in_f
            n.add_connection(_syn(G["pf_recip_inh"], exc=False),
                             src, cells[k])
        # ---- rule dress at the PF layer ----
        # heel: stance reset -> InE + InF + the E-lamination IN (g 0.5)
        heel_in = f"HEEL_{side}"
        n.add_connection(_syn(DRESS["heel_to_ine"], exc=True),
                         heel_in, f"InE_{side}")
        n.add_connection(_syn(DRESS["heel_to_inf"], exc=True),
                         heel_in, f"InF_{side}")
        n.add_connection(_syn(DRESS["heel_to_pf_in_e"], exc=True),
                         heel_in, pf_in_e)
        # toe: dorsiflexion inhibition ONLY (5 -> TOEDF -> inh 2.749)
        toe_in, toedf = f"TOE_{side}", f"TOEDF_{side}"
        self._add(toedf, TAU["preset"], n)
        n.add_connection(_syn(DRESS["toe_to_indf"], exc=True), toe_in, toedf)
        df_chan = tab["df_channel"]
        if df_chan is not None:
            n.add_connection(_syn(DRESS["df_inh"], exc=False), toedf,
                             f"PF_S{df_chan + 1}_{side}")
        # Ib load: the side LOAD IN collects stance Ib and drives
        # RG-E/InE/the E channels - gated on the EXISTING G["ib_rge"]
        # key (default 0 = edges absent). Measured in the full runner
        # air run (gate c2): force-proportional tonic Ib from the
        # standing muscle tone E-latches the RG when this reflex is
        # unconditional - the stock stack gates the same pathway on
        # ib_rge, and the runner scales the LOAD_c port with it too.
        lbin = f"LBIN_{side}"
        if G["ib_rge"] > 0.0:
            n.add_connection(_syn(DRESS["lbin_to_rge"], exc=True), lbin,
                             f"RG_E_{side}")
            n.add_connection(_syn(DRESS["lbin_to_ine"], exc=True), lbin,
                             f"InE_{side}")
            for k in range(N_SYN):
                if fam[k] == "E":
                    n.add_connection(_syn(DRESS["lbin_to_pf_e"], exc=True),
                                     lbin, cells[k])

    def _build_commissurals_side(self, n: Network, side: str):
        """Shevtsova V-class cells + their IPSILATERAL drive edges; the
        crossed edges are wired in __post_init__ after both sides exist."""
        if not self.interleg:
            return
        v2a, v0v, v0d, v3e = (f"V2a_{side}", f"V0V_{side}",
                              f"V0D_{side}", f"V3E_{side}")
        for name in (v2a, v0v, v0d, v3e):
            self._add(name, TAU["rg"], n)
        rg_e, rg_f = f"RG_E_{side}", f"RG_F_{side}"
        n.add_connection(_syn(DRESS["rgf_to_v2a"], exc=True), rg_f, v2a)
        n.add_connection(_syn(DRESS["rgf_to_v0d"], exc=True), rg_f, v0d)
        n.add_connection(_syn(DRESS["rge_to_v3"], exc=True), rg_e, v3e)
        n.add_connection(_syn(DRESS["v2a_to_v0v"], exc=True), v2a, v0v)

    def _add_muscle_neurons(self, n: Network, act: str, mi: MuscleInfo):
        # NOTE: the four motif INs are created HERE (creation pass) so the
        # antagonist cross-references in _wire_muscle_s6 (which runs in a
        # second pass over all muscles) always find their targets - the
        # lazy-create pattern of the stock builder cannot be used when the
        # motif needs BOTH directions of the IaIN/IBIN mutual inhibition.
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")
        self.mn_names[act] = mn
        self.aff_names[act] = {"Ia": ia, "II": ii, "Ib": ib}
        for name, tau in ((mn, TAU["mn"] * (1.0 + 0.5 * mi.biarticular)),
                          (ia, TAU["afferent"]), (ii, TAU["afferent"]),
                          (ib, 2.0 * TAU["afferent"]),
                          (f"IaIN_{act}", TAU["afferent"]),
                          (f"IIX_{act}", TAU["afferent"]),
                          (f"IIIN_{act}", TAU["afferent"]),
                          (f"IBIN_{act}", TAU["afferent"])):
            self._add(name, tau, n)
        n.add_input(mn)
        self.inputs.append("POST_" + act)
        for port, name in (("Ia", ia), ("II", ii), ("Ib", ib)):
            n.add_input(name)
            self.inputs.append(port + "_" + act)

    def _wire_muscle_s6(self, n: Network, act: str, mi: MuscleInfo,
                        tab: dict | None):
        mn, ia, ii, ib = (f"MN_{act}", f"Ia_{act}", f"II_{act}", f"Ib_{act}")
        d = DRESS
        covered = tab is not None and act in tab["names"]

        # ---- synergy PF -> MN (Eq-18 conductances; W > 0 only) ----
        if covered:
            row = tab["names"].index(act)
            k_vals = tab["W"][row, :]
            g_vals, valid = _eq18_conductances(k_vals)
            gs = []
            for k in range(N_SYN):
                if k_vals[k] > 0.0 and valid[k]:
                    g = float(g_vals[k])
                    gs.append(g)
                    n.add_connection(_syn(g, exc=True),
                                     f"PF_S{k + 1}_{mi.side}", mn)
            if gs:
                self.syn_g_stats.setdefault(mi.side, []).extend(gs)

        # ---- posture tonic drive (stock table + overrides) ----
        w_post = SNSpinal._group_weight(mi, SNSpinal.W_POSTURE,
                                        SNSpinal.POSTURE_OVERRIDE)
        if w_post > 0.0:
            n.add_connection(
                _syn(G["posture_to_mn"] * w_post, exc=True), "POSTURE", mn)

        # ---- Shinohara autogenic motif (rule-file gains; the four INs
        # were created in _add_muscle_neurons) ----
        n.add_connection(_syn(d["ia_to_mn"], exc=True), ia, mn)
        iain, iix, iiin, ibin = (f"IaIN_{act}", f"IIX_{act}",
                                 f"IIIN_{act}", f"IBIN_{act}")
        n.add_connection(_syn(d["ia_to_iain"], exc=True), ia, iain)
        n.add_connection(_syn(d["ii_to_iix"], exc=True), ii, iix)
        n.add_connection(_syn(d["iix_to_mn"], exc=True), iix, mn)
        n.add_connection(_syn(d["ii_to_iiin"], exc=True), ii, iiin)
        n.add_connection(_syn(d["ib_to_ibin"], exc=True), ib, ibin)
        n.add_connection(_syn(d["ibin_to_mn"], exc=False), ibin, mn)

        # antagonist projections (same side, antagonist primary groups)
        for ant in SNSpinal.ANTAGONIST.get(mi.groups[0], ()):
            for act2, mi2 in self.muscles.items():
                if mi2.side == mi.side and mi2.groups[0] == ant:
                    n.add_connection(_syn(d["iain_to_ant_mn"], exc=False),
                                     iain, f"MN_{act2}")
                    n.add_connection(_syn(d["iiin_to_ant_mn"], exc=False),
                                     iiin, f"MN_{act2}")
                    # IaIN <-> antagonist IaIN and IBIN <-> IBIN (each
                    # ordered pair wired once via this per-source loop)
                    n.add_connection(_syn(d["iain_mut"], exc=False),
                                     iain, f"IaIN_{act2}")
                    n.add_connection(_syn(d["ibin_mut"], exc=False),
                                     ibin, f"IBIN_{act2}")

        # ---- stance load sharing (extensor groups; rule-file gains) ----
        if mi.groups[0] in SNSpinal.EXTENSOR_STANCE_GROUPS:
            gname = f"IBEXC_{mi.groups[0]}_{mi.side}"
            if gname not in self.idx:
                self._add(gname, TAU["ib_exc"], n)
                self.ib_exc_groups[mi.side] = \
                    self.ib_exc_groups.get(mi.side, ()) + (mi.groups[0],)
            n.add_connection(_syn(d["ib_to_ibexc"], exc=True), ib, gname)
            n.add_connection(_syn(d["ibexc_to_mn"], exc=True), gname, mn)
            # the side LOAD IN collects stance Ib (LBIN -> RG-E/InE/PF-E
            # wired once in _build_pf_s6) - gated on G["ib_rge"] > 0
            # (default 0: force-proportional tonic Ib E-latches the RG
            # in air; see _build_pf_s6 note)
            if G["ib_rge"] > 0.0:
                n.add_connection(_syn(d["ib_to_lbin"], exc=True), ib,
                                 f"LBIN_{mi.side}")


def build(model_actuators: list[str], dt: float = SNSpinal.DT,
          interleg: bool = True) -> Syn6Network:
    """Variant build: same contract as build_network.build."""
    muscles: dict[str, MuscleInfo] = {}
    for act in model_actuators:
        mi = classify(act)
        if mi is None:
            raise ValueError(f"actuator {act!r} not in muscle_map")
        muscles[act] = mi
    sides = tuple(sorted({mi.side for mi in muscles.values()}))
    net = Syn6Network(muscles=muscles, sides=sides, interleg=interleg)
    net.compile(dt=dt)
    return net
