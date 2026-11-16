"""Export and reopen a display-only pilot route on the robotbody model."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import opensim as osim


def export(record_file: Path, geometry_dir: Path) -> Path:
    record = json.loads(record_file.read_text(encoding="utf-8"))
    reference = Path(record["reference"])
    output = record_file.with_name("tibialis_posterior_preview.osim")
    osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry_dir))
    model = osim.Model(str(reference))
    model.setName(f"tibialis_posterior_{record['branch']}_length_pilot")
    muscles = model.updMuscles()
    for i in range(muscles.getSize()):
        path = muscles.get(i).updGeometryPath()
        path.upd_Appearance().set_color(osim.Vec3(0.55, 0.55, 0.55))
    overlay = osim.PathActuator()
    overlay.setName("tibialis_posterior_blue_preview")
    overlay.set_appliesForce(False)
    overlay.updGeometryPath().upd_Appearance().set_color(osim.Vec3(0.03, 0.3, 0.85))
    for i, (body, point) in enumerate(zip(record["bodies"], record["points_m"], strict=True)):
        overlay.addNewPathPoint(f"route_{i+1}", model.getBodySet().get(body), osim.Vec3(*point))
    model.addForce(overlay)
    model.finalizeConnections()
    model.printToXML(str(output))

    # Resolve standard mesh files on easteregg2 so the preview can be
    # opened from the run folder without copying or editing the reference.
    xml = ET.parse(output)
    for mesh in xml.findall(".//mesh_file"):
        if mesh.text and (geometry_dir / mesh.text).is_file():
            mesh.text = (geometry_dir / mesh.text).as_posix()
    xml.write(output, encoding="utf-8", xml_declaration=True)

    reopened = osim.Model(str(output))
    state = reopened.initSystem()
    saved = osim.PathActuator.safeDownCast(reopened.getForceSet().get(overlay.getName()))
    points = saved.getGeometryPath().getPathPointSet()
    assert points.getSize() == len(record["bodies"])
    assert not saved.get_appliesForce()
    for i, (body, xyz) in enumerate(zip(record["bodies"], record["points_m"], strict=True)):
        point = points.get(i)
        assert point.getParentFrame().getName() == body
        loc = point.getLocation(state)
        assert max(abs(loc.get(j) - xyz[j]) for j in range(3)) < 1e-10
    print(f"Reopened {output} and verified {points.getSize()} frames and points.")
    return output


def export_edl_pair(record_file: Path, geometry_dir: Path) -> Path:
    record = json.loads(record_file.read_text(encoding="utf-8"))
    output = record_file.with_name("edl_independent_pair_preview.osim")
    routes = [("edl_1", record["routeParts"][0], (0, .19, .87)),
              ("edl_2", record["routeParts"][1], (.22, .46, .91)),
              ("ehl_20mm", record["innerRoute"], (.40, .40, .40))]
    return export_routes(record["model"], output, geometry_dir, routes)


def export_ankle20(folder: Path, geometry_dir: Path, branch: str) -> Path:
    records = [json.loads(p.read_text(encoding="utf-8"))
               for p in sorted(folder.glob(f"*_{branch}.json"))]
    assert len(records) == 10
    routes = []
    for r in records:
        if branch=="compliant" and r["group"]=="ext_dig":
            pair_file=folder/"edl_fixed_pair_separate_guides"/"ext_dig_compliant.json"
            pair=json.loads(pair_file.read_text(encoding="utf-8"))
            routes.extend((f"ext_dig_{i}_20mm_preview",part,(0,.19,.87))
                          for i,part in enumerate(pair["routeParts"],start=1))
        else:
            routes.append((r["group"]+"_20mm_preview",r,(0,.19,.87)))
    return export_routes(records[0]["model"],folder/f"ankle20_{branch}_preview.osim",
                         geometry_dir,routes)


def export_routes(reference, output, geometry_dir, routes):
    osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry_dir))
    model = osim.Model(reference)
    model.setName(output.stem)
    muscles = model.updMuscles()
    for i in range(muscles.getSize()):
        muscles.get(i).updGeometryPath().upd_Appearance().set_color(osim.Vec3(.55, .55, .55))
    for name, route, color in routes:
        overlay = osim.PathActuator()
        overlay.setName(name)
        overlay.set_appliesForce(False)
        overlay.updGeometryPath().upd_Appearance().set_color(osim.Vec3(*color))
        for i, (body, point) in enumerate(zip(route["bodies"], route["points_m"], strict=True)):
            overlay.addNewPathPoint(f"route_{i+1}", model.getBodySet().get(body), osim.Vec3(*point))
        model.addForce(overlay)
    model.finalizeConnections()
    model.printToXML(str(output))
    xml = ET.parse(output)
    for mesh in xml.findall(".//mesh_file"):
        if mesh.text and (geometry_dir / mesh.text).is_file():
            mesh.text = (geometry_dir / mesh.text).as_posix()
    xml.write(output, encoding="utf-8", xml_declaration=True)
    reopened = osim.Model(str(output))
    state = reopened.initSystem()
    n = 0
    for name, route, _ in routes:
        path = osim.PathActuator.safeDownCast(reopened.getForceSet().get(name))
        points = path.getGeometryPath().getPathPointSet()
        assert not path.get_appliesForce() and points.getSize() == len(route["bodies"])
        for i, (body, xyz) in enumerate(zip(route["bodies"], route["points_m"], strict=True)):
            point = points.get(i)
            assert point.getParentFrame().getName() == body
            assert max(abs(point.getLocation(state).get(j) - xyz[j]) for j in range(3)) < 1e-10
            n += 1
    print(f"Reopened {output} and verified {n} points on their native frames.")
    return output


if __name__ == "__main__":
    geometry = Path(sys.argv[1])
    for path in sys.argv[2:]:
        record_file = Path(path)
        if record_file.is_dir():
            for branch in ("rigid","compliant"):
                export_ankle20(record_file, geometry, branch)
        elif "routeParts" in json.loads(record_file.read_text(encoding="utf-8")):
            export_edl_pair(record_file, geometry)
        else:
            export(record_file, geometry)
