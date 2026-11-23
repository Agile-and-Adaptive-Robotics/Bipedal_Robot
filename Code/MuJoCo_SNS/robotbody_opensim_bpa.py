"""Forward integration on the mirrored Gait2392 robotbody's native OpenSim geometry."""
import xml.etree.ElementTree as ET

import numpy as np
import opensim as osim


def prepare(source, output, assignments, joints, survey_sizing=True):
    """Preserve original GeometryPaths, conditional points, joints, and inertias."""
    tree = ET.parse(source)
    model = tree.find("Model")
    forces = model.find("./ForceSet/objects")
    if len(forces) != 92 or any(f.tag != "Thelen2003Muscle" for f in forces):
        raise ValueError("Unexpected source force elements in the mirrored robotbody.")
    for old in list(forces):
        name = old.get("name")
        if name in assignments:
            actuator = ET.SubElement(forces, "PathActuator", name=name)
            actuator.append(old.find("GeometryPath"))
            ET.SubElement(actuator, "optimal_force").text = "1"
        forces.remove(old)
    if {f.get("name") for f in forces} != set(assignments):
        raise ValueError("A registered robot route is absent from the native OpenSim model.")
    pelvis = []
    for coordinate in model.findall(".//Coordinate"):
        name = coordinate.get("name")
        if name.startswith("pelvis_"):
            coordinate.find("locked").text = "true"
            pelvis.append(name)
        elif name in ("knee_angle_r", "knee_angle_l"):
            coordinate.find("default_value").text = str(np.radians(-30.))
        if name in joints:
            # In the source, range is metadata. Passive end stops let every joint
            # move inside its original bounds when the body is driven dynamically.
            lower, upper = np.degrees(np.fromstring(coordinate.findtext("range"), sep=" "))
            stop = ET.SubElement(forces, "CoordinateLimitForce", name=name + "_range_stop")
            stiffness, damping, transition = (400., .25, .1) if survey_sizing else (80., .06, .3)
            properties = dict(coordinate=name, upper_limit=upper, lower_limit=lower,
                              upper_stiffness=stiffness, lower_stiffness=stiffness,
                              damping=damping, transition=transition)
            for prop, value in properties.items():
                ET.SubElement(stop, prop).text = str(value)
    if len(pelvis) != 6:
        raise ValueError("Expected six pelvis coordinates to fix in the suspended experiment.")
    # The source has no foot contact. Gravity and passive range stops remain
    # distinct from the 54 routed pneumatic actuators.
    model.set("name", "gait2392_robotbody_bpa")
    tree.write(output, encoding="utf-8", xml_declaration=True)
    native = osim.Model(str(output))
    state = native.initSystem()
    acts = actuators(native)
    if (native.getForceSet().getSize() != 54 + len(joints) or len(acts) != 54
            or native.getControllerSet().getSize()):
        raise ValueError("Unexpected auxiliary forces or controllers in the native OpenSim model.")
    return native, [a.getName() for a in acts], np.array([a.getLength(state) for a in acts])


def actuators(model):
    routes = model.getActuators()
    return [osim.PathActuator.safeDownCast(routes.get(i)) for i in range(routes.getSize())]


def simulate(model, experiment, steps, dt, joints):
    state = model.initSystem()
    model.realizeVelocity(state)
    routes = actuators(model)
    coordinates = [model.getCoordinateSet().get(name) for name in joints]
    knee = model.getCoordinateSet().get("knee_angle_r")
    moment_arms = {a.getName(): float(a.computeMomentArm(state, knee))
                   for a in routes if a.getName() in ("vas_med_r", "bifemsh_r")}
    if not (moment_arms["vas_med_r"] > 0 > moment_arms["bifemsh_r"]):
        raise RuntimeError("Native OpenSim right-knee path signs differ from the expected anatomy.")
    for actuator in routes:
        actuator.overrideActuation(state, True)
    angles = np.empty((steps + 1, len(coordinates)))
    angles[0] = [coordinate.getValue(state) for coordinate in coordinates]
    max_error = 0.
    # Apply the actual positive tensile force in N to each GeometryPath. Starting
    # each integration from the updated state invalidates Simbody's previous forces.
    for i in range(steps):
        model.realizeVelocity(state)
        lengths = np.array([a.getLength(state) for a in routes])
        velocities = np.array([a.getLengtheningSpeed(state) for a in routes])
        forces = experiment.step(i, lengths, velocities)
        for actuator, force in zip(routes, forces):
            actuator.setOverrideActuation(state, float(force))
        manager = osim.Manager(model)
        manager.setIntegratorAccuracy(1e-7)
        manager.setWriteToStorage(False)
        manager.setPerformAnalyses(False)
        manager.initialize(state)
        target_time = (i + 1) * dt
        state = osim.State(manager.integrate(target_time))
        model.realizeDynamics(state)
        native_force = np.array([a.getActuation(state) for a in routes])
        max_error = max(max_error, float(np.max(abs(native_force - forces))))
        angles[i + 1] = [coordinate.getValue(state) for coordinate in coordinates]
        if not np.isclose(state.getTime(), target_time, rtol=0, atol=1e-12) or not np.isfinite(angles[i + 1]).all():
            raise RuntimeError(f"Native OpenSim integration failed at sample {i + 1}.")
    if max_error > 1e-9:
        raise RuntimeError("Native OpenSim route actuation disagrees with the BPA force.")
    experiment.trace.update(time_s=np.arange(steps + 1) * dt, joint_angle_rad=angles,
                            joint_names=np.array(joints))
    return experiment.trace, {"max_tensile_force_error_N": max_error, "max_ground_contacts": 0,
                              "initial_native_right_knee_moment_arm_m": moment_arms}
