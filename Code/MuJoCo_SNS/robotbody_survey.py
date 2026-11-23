"""Transfer the Appendix C 27-site BPA sizing to the mirrored Gait2392 robotbody."""
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DESIGNS = HERE.parent / "Matlab/Mesh_Optimization/Results/KompareXXVII/route_configurations.csv"
REFERENCE = HERE / "robotbody_survey_reference.csv"
ALIASES = {"vas_int": "vas_med", "glut_max": "glut_max3", "add_mag": "add_mag3"}


def load_design(variant):
    with DESIGNS.open(newline="") as source:
        rows = [r for r in csv.DictReader(source) if r["side"] == variant]
    with REFERENCE.open(newline="") as source:
        references = {r["site"]: r for r in csv.DictReader(source) if r["design"] == variant}
    if len(rows) != 27 or len(references) != 27:
        raise ValueError("The sizing and its source pose must cover 27 unique sites per side.")
    designs = {}
    for row in rows:
        source_name = row["site"][:-2]
        route_name = ALIASES.get(source_name, source_name)
        designs[route_name] = (row, references[row["site"]])
    counts = sum(int(r["nBPA"]) for r in rows)
    if ((variant == "baseline" and counts != 51) or (variant == "selected" and counts != 69)
            or len(designs) != 27):
        raise ValueError("Unexpected Appendix C actuator sizing or site aliases.")
    return designs


def configuration(names, lengths, assignments, variant, mount, reference_contraction, route_contractions=None):
    """Fix the documented sleeve, free length, count, diameter, fittings, and tendon.

    The optional transfer mount introduces a separately recorded extension
    or fixed path span. Its operating-point choice never changes the source hardware.
    The literal and survey-datum alternatives expose direct transfer incompatibilities.
    """
    designs = load_design(variant)
    result = {"command_mode": "direct_pressure", "design": variant, "mount": mount,
              "reference_contraction": reference_contraction, "actuators": {}}
    for i, name in enumerate(names):
        family = name[:-2]
        row, reference = designs[family]
        if int(row["actuator_id"]) != assignments[name]["actuator_id_per_side"]:
            raise ValueError(f"Appendix C/CAD site mapping differs for {name}.")
        diameter, count, gain = (int(row[key]) for key in ("diameter_saved_mm", "nBPA", "G_survey"))
        if int(row["tackleLineParts_c1"]) != gain or int(row["tackleLineParts_c2"]) != 1:
            raise ValueError(f"Unsupported distal or compound tackle at {name}.")
        rest, free, tendon, fitting, xi0 = (float(row[key]) for key in
            ("resting_length_m", "kmax", "tendon_length_m", "fitting_length_m", "Xi0_m"))
        fixed = tendon + 2 * fitting + xi0
        bracket = float(reference["bracket_single_m_N"])
        cable_k = float(reference["shared_tendon_N_m"])
        opposite = family + ("_l" if name.endswith("_r") else "_r")
        paired_output = float((lengths[i] + lengths[names.index(opposite)]) / 2)
        source_output = float(reference["output_route_m"])
        source_input = float(reference["input_route_m"])
        start_contraction = (route_contractions or {}).get(family, reference_contraction)
        if mount == "literal":
            equivalent_input = source_input + (paired_output - source_output) / gain
        elif mount == "survey-datum":
            equivalent_input = source_input
        else:
            equivalent_input = fixed + rest - start_contraction * (rest - free)
        # At a selected primary reverse tackle, output travel = G * sleeve
        # travel and output tension = bundle tension / G. Its datum determines
        # the initial amount of line taken up by the transmission.
        extra_output_span = gain * (equivalent_input - source_input) + (source_output - paired_output)
        if mount == "literal":
            extra_output_span = 0.
        extension = max(0., -extra_output_span)
        # A fixed span upstream of the tackle is expressed on the sleeve
        # side. One metre there takes up G metres of output line.
        extra_proximal_span = max(0., extra_output_span) / gain
        # Any added cable is on the one shared output line. Apply its
        # compliance in series with that cable from the source survey.
        output_cable_k = cable_k * tendon / (tendon + extension)
        # Each parallel sleeve has its own projected bracket compliance.
        # Their bundle drives a single output cable through a lossless tackle.
        k_effective = 1. / (bracket / count + 1. / (gain * gain * output_cable_k))
        force_source = ("manufacturer surface/capacity" if diameter == 40
                        else "Chapter 3 pressure-contraction and F620 functions")
        result["actuators"][name] = dict(
            source_site=row["site"], design_id=int(row["actuator_id"]), diameter_mm=diameter,
            resting_length_m=rest, free_length_at_620kPa_m=free,
            maximum_contraction_fraction=(rest-free)/rest, nonmuscle_length_m=fixed,
            source_tendon_length_m=tendon, source_fitting_each_m=fitting, source_Xi0_m=xi0,
            parallel_count=count, transmission_gain=gain, source_Xi1_N_m=float(row["Xi1_N_per_m"]),
            source_Xi2_N_m=float(row["Xi2_N_per_m"]), single_bracket_compliance_m_N=bracket,
            source_output_cable_stiffness_N_m=cable_k, output_cable_stiffness_N_m=output_cable_k,
            series_stiffness_N_m=k_effective, source_reference_output_m=source_output,
            source_reference_input_m=source_input, output_reference_m=paired_output,
            bpa_side_reference_m=equivalent_input, extra_fixed_proximal_span_m=extra_proximal_span,
            extra_output_cable_m=extension, pressure_max_kpa=float(row["pressure_kPa"]),
            transfer_reference_contraction=start_contraction,
            tau_fill_s=.040, tau_vent_s=.060, delay_s=0.,
            provenance=f"Appendix C {variant} sizing, {row['site']}; {force_source}. "
                       "Saved design free length; valve timing assumed. See configuration mount and source MAT.",
            source_mat=row["mat_file"].replace("\\","/").rsplit("/",1)[-1],
        )
    return result
