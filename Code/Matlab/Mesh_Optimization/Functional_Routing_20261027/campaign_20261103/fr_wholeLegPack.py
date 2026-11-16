"""Search coordinated sleeve lanes and verify the chosen whole-leg layouts."""
import csv
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pad(option):
    poses = option["parts"]
    a = np.zeros((len(poses), max(map(len, poses)), 9))
    a[:, :, 7] = -1
    for i, parts in enumerate(poses):
        a[i, :len(parts)] = parts
    return a


def distances(a, b):
    """Exact distances between finite 3-D line segments, including points."""
    u, v = a[:, 3:6] - a[:, :3], b[:, 3:6] - b[:, :3]
    w = a[:, :3] - b[:, :3]
    dot = lambda x, y: np.einsum("ij,ij->i", x, y)
    uu, uv, vv = dot(u, u), dot(u, v), dot(v, v)
    uw, vw = dot(u, w), dot(v, w)
    us, vs = np.maximum(uu, 1e-30), np.maximum(vv, 1e-30)
    det = uu * vv - uv * uv
    ds = np.maximum(det, 1e-30)
    s, t = (uv * vw - vv * uw) / ds, (uu * vw - uv * uw) / ds
    inside = (det > 1e-16 * uu * vv) & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
    d = np.full(len(a), np.inf)
    for s_edge, t_edge in ((0, vw / vs), (1, (uv + vw) / vs),
                           (-uw / us, 0), ((uv - uw) / us, 1),
                           (np.where(inside, s, 0), np.where(inside, t, 0))):
        e = w + np.asarray(np.clip(s_edge, 0, 1))[..., None] * u
        e -= np.asarray(np.clip(t_edge, 0, 1))[..., None] * v
        d = np.minimum(d, dot(e, e))
    return np.sqrt(d)


class Envelope:
    def __init__(self, option):
        self.a = pad(option)
        r = self.a[:, :, 6, None]
        self.lo = np.minimum(self.a[:, :, :3], self.a[:, :, 3:6]) - r
        self.hi = np.maximum(self.a[:, :, :3], self.a[:, :, 3:6]) + r
        self.bounds = (self.lo.min(axis=(0, 1)), self.hi.max(axis=(0, 1)))


def clearance(a, b, same=False):
    # Values farther than 18 mm are recorded as a lower bound.
    limit = .018
    if np.any(a.bounds[0] - b.bounds[1] > limit) or np.any(b.bounds[0] - a.bounds[1] > limit):
        return (18., 18., 18.), (-1, -1, -1)
    keep = (a.a[:, :, 7, None] >= 0) & (b.a[:, None, :, 7] >= 0)
    if same:
        keep &= a.a[:, :, 8, None] != b.a[:, None, :, 8]
        keep &= np.arange(a.a.shape[1])[:, None] < np.arange(b.a.shape[1])
    for j in range(3):
        keep &= (a.lo[:, :, j, None] - b.hi[:, None, :, j] <= limit)
        keep &= (b.lo[:, None, :, j] - a.hi[:, :, j, None] <= limit)
    pose, u, v = np.nonzero(keep)
    if not len(pose):
        return (18., 18., 18.), (-1, -1, -1)
    x, y = a.a[pose, u], b.a[pose, v]
    gap = 1000 * (distances(x, y) - x[:, 6] - y[:, 6])
    kind_a, kind_b = x[:, 7], y[:, 7]
    is_sleeve = (kind_a == 0) & (kind_b == 0)
    is_hardware = (kind_a == 1) | (kind_b == 1)
    mins, witnesses = [], []
    for mask in (is_sleeve, is_hardware, ~(is_sleeve | is_hardware)):
        candidates = np.flatnonzero(mask)
        if len(candidates):
            i = candidates[np.argmin(gap[candidates])]
            mins.append(min(18., float(gap[i])))
            witnesses.append(int(pose[i]) if gap[i] < 18 else -1)
        else:
            mins.append(18.)
            witnesses.append(-1)
    return tuple(mins), tuple(witnesses)


def conflict_score(gaps):
    sleeve, hardware, cable = gaps
    return (70 * (sleeve < 0) + 4 * max(0, -sleeve)
            + 90 * (hardware < 0) + 4 * max(0, -hardware)
            + 15 * (cable < 0) + 2 * max(0, -cable)
            + .12 * min(max(sleeve, 0), 12))


def unpack(folder, branch):
    options = {}
    for file in sorted(folder.glob(f"*_{branch}.json")):
        variants = read(file)
        options[variants[0]["group"]] = {v["id"]: v for v in variants}
    assert len(options) == 27
    return options


def rank_local(option, baseline, internal):
    limits = [option.get(k) for k in ("bone_mm", "foot_lead_mm", "home_side_mm",
                                      "sleeve_bay_mm", "tendon_span_mm", "foot_occupancy_mm")]
    hard = sum(75 + 5 * max(0, -v) for v in limits if v is not None and v < 0)
    torque = 2.5 * max(0, baseline["capacity"] - option["capacity"])
    force = .7 * max(0, baseline["active"] - option["active"])
    motion = .15 * np.linalg.norm(option["offset"][:3])
    return hard + torque + force + motion + conflict_score(internal)


def solve(options, seed):
    groups = sorted(options)
    shapes, costs, states = [], [], []
    for name in groups:
        variants = options[name]
        env = {i: Envelope(r) for i, r in variants.items()}
        cost = {i: rank_local(r, variants[1], clearance(env[i], env[i], same=True)[0])
                for i, r in variants.items()}
        order = sorted(variants, key=cost.get)
        # Keep the current route and diverse offsets, along with good local scores.
        ids = list(dict.fromkeys([1] + order[:11]))
        shapes.append(env)
        costs.append(cost)
        states.append(ids)

    @lru_cache(maxsize=None)
    def cached(g1, o1, g2, o2):
        return conflict_score(clearance(shapes[g1][o1], shapes[g2][o2])[0])

    def pair(i, a, j, b):
        return cached(i, a, j, b) if i < j else cached(j, b, i, a)

    def total(chosen):
        return sum(costs[i][a] for i, a in enumerate(chosen)) + sum(
            pair(i, chosen[i], j, chosen[j]) for i in range(len(groups)) for j in range(i))

    rng = np.random.default_rng(seed)
    winner, value = [1] * len(groups), np.inf
    for attempt in range(5):
        selected = winner.copy()
        if attempt:
            for i in rng.choice(len(groups), 9, replace=False):
                selected[i] = int(rng.choice(states[i]))
        for _ in range(5):
            changed = False
            for i in rng.permutation(len(groups)):
                i = int(i)
                old = selected[i]
                def local(a):
                    return costs[i][a] + sum(pair(i, a, j, b) for j, b in enumerate(selected) if j != i)
                selected[i] = min(states[i], key=local)
                changed |= selected[i] != old
            if not changed:
                break
        score = total(selected)
        if score < value:
            winner, value = selected, score
        print(f"Layout restart {attempt + 1}: score {score:.1f}, cached pairs {cached.cache_info().currsize}",
              flush=True)
    print(f"Selected score {value:.1f} versus unchanged score {total([1] * len(groups)):.1f}", flush=True)
    return dict(zip(groups, winner))


def verify(options, selected, poses):
    groups = sorted(options)
    env = {g: Envelope(options[g][selected[g]]) for g in groups}
    rows, per_group = [], {}
    for i, g in enumerate(groups):
        r = options[g][selected[g]]
        self_gaps, _ = clearance(env[g], env[g], same=True)
        per_group[g] = dict(group=g, option=selected[g], actuators=sum(x["count"] for x in r["lines"]),
                            bone_mm=r["bone_mm"], foot_lead_mm=r["foot_lead_mm"], home_side_mm=r["home_side_mm"],
                            sleeve_bay_mm=r.get("sleeve_bay_mm"), tendon_span_mm=r.get("tendon_span_mm"),
                            foot_occupancy_mm=r.get("foot_occupancy_mm"),
                            sleeve_gap_mm=self_gaps[0], hardware_gap_mm=self_gaps[1], cable_gap_mm=self_gaps[2],
                            signed_capacity_pct=r["capacity"], target_met_pct=r["target_met"],
                            active_force_pct=r["active"], common_active_pct=r["common_active"], peak_tendon_N=r["peak_N"])
        for h in groups[i + 1:]:
            gaps, hits = clearance(env[g], env[h])
            rows.append(dict(group_a=g, group_b=h, sleeve_gap_mm=gaps[0], hardware_gap_mm=gaps[1],
                             cable_gap_mm=gaps[2], sleeve_pose=hits[0], hardware_pose=hits[1], cable_pose=hits[2]))
            for name in ("sleeve", "hardware", "cable"):
                per_group[g][name + "_gap_mm"] = min(per_group[g][name + "_gap_mm"], gaps[("sleeve", "hardware", "cable").index(name)])
    for r in rows:
        for key in ("sleeve_gap_mm", "hardware_gap_mm", "cable_gap_mm"):
            per_group[r["group_b"]][key] = min(per_group[r["group_b"]][key], r[key])
    for g, r in per_group.items():
        keys = ("bone_mm", "foot_lead_mm", "home_side_mm", "sleeve_bay_mm",
                "tendon_span_mm", "foot_occupancy_mm", "sleeve_gap_mm", "hardware_gap_mm", "cable_gap_mm")
        failures = [key for key in keys if r[key] is not None and r[key] < -1e-6]
        r["clearance_pass"] = not failures
        if r["active_force_pct"] < 100 or r["common_active_pct"] < 100:
            failures.append("no_force_poses")
        r["failed_checks"] = ", ".join(failures)
    return list(per_group.values()), rows


def main(folder, phase):
    branches = ("rigid", "compliant")
    if phase == "select":
        choices = {}
        for i, b in enumerate(branches):
            choices.update({f"{g}_{b}": int(index)
                            for g, index in solve(unpack(folder / "coarse", b), 52 + i).items()})
        write(folder / "selection.json", choices)
        return
    poses = read(folder / "dense" / "poses.json")
    choices = read(folder / "selection.json")
    summaries, all_groups = [], {}
    from fr_exportPilotOpenSim import export_routes
    for b in branches:
        options = unpack(folder / "dense", b)
        selected = {g: choices[f"{g}_{b}"] for g in options}
        for label, pick in (("before", dict.fromkeys(options, 1)), ("packed", selected)):
            groups, pairs = verify(options, pick, poses)
            table(folder / f"{b}_{label}_groups.csv", groups)
            table(folder / f"{b}_{label}_pairs.csv", pairs)
            summary = dict(branch=b, layout=label, poses=len(poses["angles_deg"]),
                           clearance_passes=sum(r["clearance_pass"] for r in groups),
                           sleeve_clear_groups=sum(r["sleeve_gap_mm"] >= 0 for r in groups),
                           hardware_clear_groups=sum(r["hardware_gap_mm"] >= 0 for r in groups),
                           sleeve_pairs=sum(r["sleeve_gap_mm"] < 0 for r in pairs),
                           hardware_pairs=sum(r["hardware_gap_mm"] < 0 for r in pairs),
                           cable_pairs=sum(r["cable_gap_mm"] < 0 for r in pairs),
                           **{key: float(np.mean([r[key] for r in groups]))
                              for key in ("signed_capacity_pct", "target_met_pct", "active_force_pct", "common_active_pct")})
            summaries.append(summary)
            all_groups[b, label] = groups
            print(summary, flush=True)
        routes = [(f"{g}_{r['channel']}", r, (0, .22, .87)) for g, opt in options.items()
                  for r in opt[selected[g]]["lines"]]
        export_routes(poses["reference"], folder / f"whole_leg_{b}_packed.osim",
                      Path("C:/OpenSim 4.6/Geometry"), routes)
    table(folder / "packaging_summary.csv", summaries)
    report(folder, summaries, all_groups)


def report(folder, summaries, all_groups):
    lines = ["# Whole-leg packaging verification", "",
             "The search keeps the 27 right-side groups and their checkpoint actuator sizes and counts. "
             "The ankle starts from the previously coordinated separate-line geometry. Each search "
             "alternative moves native-frame mounts, staggers fittings, or changes bundle spread. "
             "Foot attachments stay fixed. Torque and inflation are solved again after changes.", "",
             "| Branch | Groups clearing modeled checks | Sleeve-conflict pairs | Hardware-conflict pairs | "
             "Signed target poses met | Active force on primary sweep | Signed capacity retained |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for branch in ("rigid", "compliant"):
        b, a = [s for s in summaries if s["branch"] == branch]
        retention = 100 * a["signed_capacity_pct"] / max(b["signed_capacity_pct"], 1e-9)
        lines.append(f"| {branch.title()} | {b['clearance_passes']} → {a['clearance_passes']} of 27 | "
                     f"{b['sleeve_pairs']} → {a['sleeve_pairs']} | {b['hardware_pairs']} → {a['hardware_pairs']} | "
                     f"{b['target_met_pct']:.1f}% → {a['target_met_pct']:.1f}% | "
                     f"{b['active_force_pct']:.1f}% → {a['active_force_pct']:.1f}% | {retention:.1f}% |")
    lines += ["", f"Both layouts were replayed on {summaries[0]['poses']} common poses, with knee motion "
              "bounded to −120° through 0°. A primary sweep also compares signed torque and force on "
              "31 points per group. Signed capacity clips each point at its human target before "
              "averaging equally across groups. Retention compares packed and starting capacity.", "",
              "The bone check uses the verified femur, tibia, and talus fields. The pelvis, trunk, "
              "calcaneus, and toes lack a verified signed surface field in this run. The CAD registration "
              "report gives station errors of 87.0 mm at the knee and 88.8 mm at the ankle. The rigid "
              "hardware envelopes are fitting, tackle, and stop proxies. OpenSim previews show nominal "
              "centerlines, with physical sleeves and cable envelopes in the packaging graphics. "
              "The packaged routes retain the general 27-group topology, including the combined vasti "
              "candidate and general gastrocnemius routes from the main checkpoints.", ""]
    for branch in ("rigid", "compliant"):
        groups = all_groups[branch, "packed"]
        clear = [r["group"] for r in groups if r["clearance_pass"]]
        sleeve_clear = [r["group"] for r in groups if r["sleeve_gap_mm"] >= 0]
        lines += [f"## {branch.title()}", "",
                  "Groups clearing the modeled checks: " + (", ".join(clear) or "none") + ".", "",
                  "Groups with no inflated sleeve overlap: " + (", ".join(sleeve_clear) or "none") + ".", "",
                  "| Group | Min. sleeve gap | Min. hardware gap | Torque target met | "
                  "Active force | Failed checks |",
                  "|---|---:|---:|---:|---:|---|"]
        for r in groups:
            lines.append(f"| {r['group']} | {r['sleeve_gap_mm']:+.1f} mm | "
                         f"{r['hardware_gap_mm']:+.1f} mm | {r['target_met_pct']:.1f}% | "
                         f"{r['active_force_pct']:.1f}% | {r['failed_checks'] or '—'} |")
        lines += ["", "Positive gaps reported at 18 mm are lower bounds. Pair CSV pose indices "
                  "are zero-based into `dense/poses.json`.", ""]
    (folder / "WHOLE_LEG_PACKAGING_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2])
