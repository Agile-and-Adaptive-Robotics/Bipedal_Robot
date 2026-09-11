"""Functional map of the Gait2392 muscle set (92 actuators in the MuJoCo MJCF).

Each actuator name (minus the _r/_l suffix) is classified into functional
groups used by the spinal network: reflex routing (which afferents gate where),
pattern-formation (PF) -> motoneuron initial weights, and posture drive.

Group membership is a list; the first entry is the PRIMARY group (used for
phase-gated drive), later entries are secondary (half weight). Muscles that
cross two joints (biarticular) get flagged automatically.

Sources for the group assignments: the Gait2392 model documentation
(Delp et al. 1990 / gait2392_thelen2003muscle.osim) and standard clinical
functional grouping; see spinal/DESIGN.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field

GROUPS = ("hip_ext", "hip_flex", "hip_abd", "hip_add",
          "knee_ext", "knee_flex", "ankle_pf", "ankle_df",
          "trunk_ext", "trunk_flex")

# name (no side suffix) -> functional groups, primary first
_GROUPS_BY_NAME: dict[str, tuple[str, ...]] = {
    # --- hip abductors (stance lateral stability) ---
    "glut_med1": ("hip_abd",),
    "glut_med2": ("hip_abd",),
    "glut_med3": ("hip_abd",),
    "glut_min1": ("hip_abd",),
    "glut_min2": ("hip_abd",),
    "glut_min3": ("hip_abd",),
    # --- hamstrings (bi: hip ext + knee flex) ---
    "semimem": ("knee_flex", "hip_ext"),
    "semiten": ("knee_flex", "hip_ext"),
    "bifemlh": ("knee_flex", "hip_ext"),
    "bifemsh": ("knee_flex",),
    # --- hip flexors ---
    "sar": ("hip_flex",),                      # biarticular in anatomy; its
                                               # knee-flex ride-along on the
                                               # strong F1 drive over-drove it
                                               # (peak 0.84, abduction component
                                               # splaying the swing leg)
    "iliacus": ("hip_flex",),
    "psoas": ("hip_flex",),
    "rect_fem": ("knee_ext", "hip_flex"),      # bi
    "tfl": ("hip_flex",),                      # fascia lata; its hip_abd
                                               # ride-along on the F1 flexor
                                               # drive abducted the swing leg
                                               # every cycle (2026-09-11)
    # --- hip adductors ---
    "add_long": ("hip_add",),
    "add_brev": ("hip_add",),
    "add_mag1": ("hip_add",),
    "add_mag2": ("hip_add",),
    "add_mag3": ("hip_ext", "hip_add"),        # extensor (vertical) fibers
    "pect": ("hip_add",),
    "grac": ("hip_add",),                      # biarticular in anatomy, but
                                               # its knee-flex ride-along on
                                               # the strong F1 drive splayed
                                               # the legs +-20 deg in the
                                               # frontal plane (2026-09-11)
    # --- hip extensors (glutes + short rotators) ---
    "glut_max1": ("hip_ext",),
    "glut_max2": ("hip_ext",),
    "glut_max3": ("hip_ext",),
    "quad_fem": ("hip_ext",),
    "gem": ("hip_ext",),
    "peri": ("hip_ext",),
    # --- knee extensors ---
    "vas_med": ("knee_ext",),
    "vas_int": ("knee_ext",),
    "vas_lat": ("knee_ext",),
    # --- triceps surae + plantar flexors ---
    "med_gas": ("ankle_pf", "knee_flex"),      # bi
    "lat_gas": ("ankle_pf", "knee_flex"),      # bi
    "soleus": ("ankle_pf",),
    "tib_post": ("ankle_pf",),
    "flex_dig": ("ankle_pf",),
    "flex_hal": ("ankle_pf",),
    "per_brev": ("ankle_pf",),
    "per_long": ("ankle_pf",),
    "per_tert": ("ankle_pf",),
    # --- dorsiflexors ---
    "tib_ant": ("ankle_df",),
    "ext_dig": ("ankle_df",),
    "ext_hal": ("ankle_df",),
    # --- trunk ---
    "ercspn": ("trunk_ext",),
    "intobl": ("trunk_flex",),
    "extobl": ("trunk_flex",),
}


@dataclass
class MuscleInfo:
    actuator: str                 # full MJCF actuator name, e.g. "vas_lat_r"
    side: str                     # "r" or "l" or "x" (trunk singles stay per-side)
    base: str                     # name without side suffix
    groups: tuple[str, ...] = field(default_factory=tuple)
    biarticular: bool = False
    fmax: float = 1.0             # peak isometric force [N] from MJCF gainprm
    length0: float = 0.0          # tendon length at the keyframe pose [m]


def classify(actuator_name: str) -> MuscleInfo | None:
    """Split an actuator name into side/base and look up its groups."""
    if actuator_name.endswith(("_r", "_l")):
        side = actuator_name[-1]
        base = actuator_name[:-2]
    else:
        side, base = "x", actuator_name
    groups = _GROUPS_BY_NAME.get(base)
    if groups is None:
        return None
    return MuscleInfo(actuator=actuator_name, side=side, base=base,
                      groups=groups, biarticular=len(groups) > 1)


def group_members(muscles: dict[str, MuscleInfo], side: str, group: str):
    """All actuator names on one side belonging to a functional group."""
    return [m.actuator for m in muscles.values()
            if m.side == side and group in m.groups]


EXTENSOR_STANCE_GROUPS = ("hip_ext", "knee_ext", "ankle_pf", "hip_abd", "trunk_ext")
