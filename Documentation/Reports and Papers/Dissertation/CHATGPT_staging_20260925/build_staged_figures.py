"""Build the reviewed dissertation figures prepared on 2026-09-25.

The Xi review copies are written into this staging directory.  The approved
Steele source figure is written to ``Figures/Aim2``; ``ProofFinal`` receives a
copy as part of the dissertation synchronization step, outside this script.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
from PIL import Image
from pypdf import PdfReader, PdfWriter


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

BLUE = "#4477AA"
ORANGE = "#EE6677"
GREEN = "#228833"
PURPLE = "#AA3377"
GRAY = "#777777"
LIGHT = "#BBBBBB"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def arrow3(ax, start, end, color="black", lw=1.7, ratio=0.10, zorder=5):
    start = np.asarray(start, dtype=float)
    d = np.asarray(end, dtype=float) - start
    ax.quiver(
        *start,
        *d,
        color=color,
        linewidth=lw,
        arrow_length_ratio=ratio,
        normalize=False,
        pivot="tail",
        zorder=zorder,
    )


def setup_3d(ax, lim=1.35, elev=23, azim=-55):
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_zlim(0, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_proj_type("ortho")
    except AttributeError:
        pass
    ax.set_axis_off()


def draw_basis(
    ax,
    R=None,
    length=1.0,
    alpha=1.0,
    labels=None,
    colors=None,
    lw=1.8,
    origin=None,
):
    if R is None:
        R = np.eye(3)
    if labels is None:
        labels = (r"$\hat{e}_1$", r"$\hat{e}_2$", r"$\hat{e}_3$")
    if colors is None:
        colors = (BLUE, ORANGE, GREEN)
    if origin is None:
        origin = np.zeros(3)
    origin = np.asarray(origin, dtype=float)
    for i, (label, color) in enumerate(zip(labels, colors)):
        end = origin + length * R[:, i]
        arrow3(ax, origin, end, color=color, lw=lw, ratio=0.11)
        label_pos = origin + 1.10 * length * R[:, i]
        ax.text(*label_pos, label, color=color, fontsize=11, alpha=alpha)


def build_projection_figure():
    fig, ax = plt.subplots(figsize=(5.8, 7.0))
    ax.set_xlim(0.0, 5.8)
    ax.set_ylim(1.1, 8.1)
    ax.set_aspect("equal")
    ax.axis("off")

    origin_1 = np.array([1.45, 2.05])
    origin_2 = np.array([4.10, 6.65])
    line = origin_2 - origin_1
    line_length = np.linalg.norm(line)
    u12 = line / line_length
    u21 = -u12
    normal = np.array([-u12[1], u12[0]])

    # Oblique projections of two Cartesian frames after rotating the former
    # viewpoint -90 degrees about e3. Every axis is also text-labeled, so the
    # figure remains interpretable without color.
    bases = (
        (
            origin_1,
            (np.array([0.86, -0.18]), np.array([-0.50, -0.58]), np.array([0.00, 0.98])),
            1,
        ),
        (
            origin_2,
            (np.array([0.65, 0.46]), np.array([-0.64, 0.20]), np.array([0.00, 0.98])),
            2,
        ),
    )
    axis_colors = (BLUE, ORANGE, GREEN)
    axis_names = ("1", "2", "3")
    for origin, basis, frame_index in bases:
        for direction, color, axis_name in zip(basis, axis_colors, axis_names):
            endpoint = origin + 0.78 * direction
            ax.annotate(
                "",
                xy=endpoint,
                xytext=origin,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.9, mutation_scale=12),
                zorder=4,
            )
            label_offset = 0.12 * unit(direction)
            ax.text(
                *(endpoint + label_offset),
                rf"$\hat e_{{{axis_name}}}^{{({frame_index})}}$",
                color=color,
                fontsize=10.5,
                ha="center",
                va="center",
            )

    # Frame names are placed in clear whitespace rather than on an axis.
    ax.text(*(origin_1 + np.array([-0.42, 0.18])), r"$\{br,1\}$", fontsize=11.5, ha="right")
    ax.text(*(origin_2 + np.array([0.42, -0.65])), r"$\{br,2\}$", fontsize=11.5, ha="left")

    # Dashed force path and equal-and-opposite forces directed toward each other.
    ax.plot(
        [origin_1[0], origin_2[0]],
        [origin_1[1], origin_2[1]],
        color="0.55",
        lw=1.0,
        ls=(0, (4, 3)),
        zorder=1,
    )
    ax.text(
        *(origin_1 + 2.05 * u12 + 0.34 * normal),
        r"$\hat{\mathbf{u}}_1=-\hat{\mathbf{u}}_2$",
        color="0.35",
        fontsize=10,
        ha="center",
        va="center",
    )
    force_length = 1.62
    for origin, direction in ((origin_1, u12), (origin_2, u21)):
        start = origin + 0.10 * direction
        end = origin + force_length * direction
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2.5, mutation_scale=14),
            zorder=6,
        )
    ax.text(
        *(origin_2 - 0.84 * u12 + 0.52 * normal),
        r"$\mathbf{F}'_2=-\mathbf{F}'_1$",
        fontsize=11.5,
        ha="center",
        va="center",
    )

    eta_color = "#882255"
    delta_color = "#332288"
    tendon_color = "#CC6677"
    eta_1 = np.array([0.58, 0.94])
    eta_2 = np.array([-0.67, 0.92])
    for origin, eta_vec, label, label_offset in (
        (origin_1, eta_1, r"$\mathbf{\eta}_{br,1}$", np.array([-0.28, 0.12])),
        (origin_2, eta_2, r"$\mathbf{\eta}_{br,2}$", np.array([-0.02, 0.17])),
    ):
        endpoint = origin + eta_vec
        ax.annotate(
            "",
            xy=endpoint,
            xytext=origin,
            arrowprops=dict(arrowstyle="-|>", color=eta_color, lw=2.7, mutation_scale=14),
            zorder=7,
        )
        ax.text(*(endpoint + label_offset), label, color=eta_color, fontsize=11.5)

    # Orthogonal construction lines connect each eta endpoint to its force-path
    # projection. The scalar deltas are drawn on offset rails with double heads.
    delta_1 = float(u12 @ eta_1)
    delta_2 = float(u21 @ eta_2)
    proj_1 = origin_1 + delta_1 * u12
    proj_2 = origin_2 + delta_2 * u21
    ax.plot(
        [origin_1[0] + eta_1[0], proj_1[0]],
        [origin_1[1] + eta_1[1], proj_1[1]],
        color="0.55",
        lw=1.1,
        ls=(0, (3, 2)),
    )
    ax.plot(
        [origin_2[0] + eta_2[0], proj_2[0]],
        [origin_2[1] + eta_2[1], proj_2[1]],
        color="0.55",
        lw=1.1,
        ls=(0, (3, 2)),
    )
    delta_1_offset = -0.23 * normal
    delta_2_offset = -0.55 * normal
    for origin, projection, offset, label, label_shift in (
        (origin_1, proj_1, delta_1_offset, r"$\delta_{br,1}$", -0.20 * normal),
        (origin_2, proj_2, delta_2_offset, r"$\delta_{br,2}$", -0.25 * normal),
    ):
        start = origin + offset
        end = projection + offset
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(
                arrowstyle="<->",
                color=delta_color,
                lw=2.8,
                mutation_scale=13,
                shrinkA=0,
                shrinkB=0,
            ),
            zorder=8,
        )
        ax.text(
            *(0.5 * (start + end) + label_shift),
            label,
            color=delta_color,
            fontsize=11.5,
            ha="center",
            va="center",
        )

    # Tendon elongation continues the bracket-1 projection on the same offset rail.
    tendon_start = proj_1 + delta_1_offset
    tendon_end = tendon_start + 0.90 * u12
    ax.annotate(
        "",
        xy=tendon_end,
        xytext=tendon_start,
        arrowprops=dict(
            arrowstyle="<->",
            color=tendon_color,
            lw=2.8,
            mutation_scale=13,
            shrinkA=0,
            shrinkB=0,
        ),
        zorder=8,
    )
    ax.text(
        *(0.5 * (tendon_start + tendon_end) - 0.34 * normal),
        r"$\delta_{\mathrm{tendon}}$",
        color=tendon_color,
        fontsize=11.5,
        ha="center",
        va="center",
    )

    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.99)
    fig.savefig(HERE / "xiProjection.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(HERE / "xiProjection.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(HERE / "xiProjection.svg", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    # Carry the long description with the standalone review PDF. The matching
    # LaTeX source also includes it as an ALT TEXT comment for Overleaf insertion.
    projection_pdf = HERE / "xiProjection.pdf"
    reader = PdfReader(str(projection_pdf))
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    alt_text = (
        "Two three-dimensional bracket coordinate frames, br comma 1 and br comma 2, are separated "
        "along a dashed actuator line. Equal-and-opposite force vectors point from each "
        "bracket toward the other and are labeled F prime 2 equals negative F prime 1. "
        "Each bracket has a local deflection vector eta and a "
        "clearly offset scalar projection delta along its force direction. The bracket 2 "
        "projection is shifted to the right of its coordinate axes. Tendon elongation "
        "delta_tendon is collinear and contiguous with the bracket 1 projection."
    )
    writer.add_metadata(
        {
            "/Title": "Two-bracket deflection and force-path projection",
            "/Subject": alt_text,
            "/AltText": alt_text,
            "/Keywords": "bracket frames, force vectors, tendon elongation, accessible scientific figure",
        }
    )
    temp_pdf = HERE / "xiProjection.metadata.pdf"
    with temp_pdf.open("wb") as stream:
        writer.write(stream)
    temp_pdf.replace(projection_pdf)


def rz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def ry(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def plot_rotation_arc(ax, points, color, label, label_offset=(0, 0, 0)):
    points = np.asarray(points)
    ax.plot(points[:, 0], points[:, 1], points[:, 2], color=color, lw=1.7)
    arrow3(ax, points[-3], points[-1], color=color, lw=1.5, ratio=0.24)
    mid = points[len(points) // 2] + np.asarray(label_offset)
    ax.text(*mid, label, color=color, fontsize=11)


def build_transform_figure():
    fig, axes = plt.subplots(2, 1, figsize=(5.7, 7.8))
    for ax in axes:
        ax.set_xlim(0.0, 6.0)
        ax.set_ylim(0.0, 3.8)
        ax.set_aspect("equal")
        ax.axis("off")

    axis_gray = "#AAAAAA"
    frame_blue = BLUE
    frame_green = GREEN
    frame_red = ORANGE
    angle_color = PURPLE

    # Panel A: Rz aligns x_1 with the sagittal-plane projection of p_A or p_B.
    ax = axes[0]
    origin = np.array([1.20, 0.78])
    xs = np.array([1.00, -0.24])
    ys = np.array([0.60, 0.52])
    zs = np.array([0.00, 1.05])
    plane = np.vstack(
        (
            origin,
            origin + 2.55 * xs,
            origin + 2.55 * xs + 1.55 * ys,
            origin + 1.55 * ys,
        )
    )
    ax.add_patch(
        patches.Polygon(plane, closed=True, facecolor="#F3F3F3", edgecolor="0.80", lw=0.9, zorder=0)
    )
    projection = origin + 1.75 * xs + 0.88 * ys
    point = projection + 0.78 * zs
    x1 = unit(projection - origin)
    y1 = np.array([-x1[1], x1[0]])

    for direction, label in ((xs, r"$\hat x_s$"), (ys, r"$\hat y_s$")):
        end = origin + 1.03 * direction
        ax.annotate("", xy=end, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=axis_gray, lw=1.4, mutation_scale=12))
        ax.text(*(end + 0.10 * unit(direction)), label, color="0.45", fontsize=10)
    ax.annotate("", xy=projection, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=frame_blue, lw=2.5, mutation_scale=14))
    ax.annotate("", xy=origin + 1.12 * zs, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=frame_green, lw=2.0, mutation_scale=12))
    ax.text(*(projection - 0.20 * y1), r"$\hat x_1$", color=frame_blue, fontsize=11)
    ax.text(*(origin + 1.22 * zs), r"$\hat z_1=\hat z_s$", color=frame_green, fontsize=11, ha="center")
    ax.plot([projection[0], point[0]], [projection[1], point[1]], color="0.45", lw=1.2, ls=(0, (3, 2)))
    ax.scatter(*projection, s=26, facecolor="white", edgecolor="black", zorder=6)
    ax.scatter(*point, s=34, color="black", zorder=6)
    ax.text(*(projection + np.array([0.14, -0.18])), r"$\mathbf{p}_{A/B}^{\mathrm{sag}}$", fontsize=10.5)
    ax.text(*(point + np.array([0.12, 0.06])), r"$\mathbf{p}_{A/B}$", fontsize=10.5)
    ax.text(3.82, 0.34, "sagittal plane", color="0.42", fontsize=9.5)
    phi0 = np.arctan2(xs[1], xs[0])
    phi1 = np.arctan2(x1[1], x1[0])
    t = np.linspace(phi0, phi1, 40)
    arc = origin + 0.72 * np.c_[np.cos(t), np.sin(t)]
    ax.plot(arc[:, 0], arc[:, 1], color=angle_color, lw=1.8)
    ax.text(*(origin + 0.86 * unit(np.array([np.cos((phi0 + phi1) / 2), np.sin((phi0 + phi1) / 2)]))), r"$\theta_z$", color=angle_color, fontsize=11)
    ax.text(0.12, 3.48, "(A)", fontsize=12, weight="bold")
    ax.text(0.62, 3.48, r"First rotation: $\hat x_1$ targets the sagittal projection", fontsize=11)

    # Panel B: Ry^T rotates about current y_1=y_br until x_br targets p_A or p_B.
    ax = axes[1]
    origin = np.array([1.18, 0.70])
    x1 = np.array([1.00, 0.00])
    y1 = np.array([-0.58, 0.42])
    z1 = np.array([0.00, 1.00])
    projection = origin + 2.22 * x1
    point = projection + 1.12 * z1
    xbr = unit(point - origin)
    zbr = np.array([-xbr[1], xbr[0]])
    plane = np.vstack((origin, origin + 2.75 * x1, origin + 2.75 * x1 + 1.45 * z1, origin + 1.45 * z1))
    ax.add_patch(
        patches.Polygon(plane, closed=True, facecolor="#F3F3F3", edgecolor="0.80", lw=0.9, zorder=0)
    )
    for direction, label in ((x1, r"$\hat x_1$"), (z1, r"$\hat z_1$")):
        end = origin + 1.12 * direction
        ax.annotate("", xy=end, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=axis_gray, lw=1.4, mutation_scale=12))
        ax.text(*(end + 0.10 * unit(direction)), label, color="0.50", fontsize=10)
    ax.annotate("", xy=point, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=frame_blue, lw=2.6, mutation_scale=14))
    ax.annotate("", xy=origin + 1.08 * y1, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=frame_red, lw=2.0, mutation_scale=12))
    ax.annotate("", xy=origin + 1.10 * zbr, xytext=origin, arrowprops=dict(arrowstyle="-|>", color=frame_green, lw=2.0, mutation_scale=12))
    xbr_normal = np.array([-xbr[1], xbr[0]])
    ax.text(*(origin + 1.72 * xbr - 0.19 * xbr_normal), r"$\hat x_{br}$", color=frame_blue, fontsize=11)
    ax.text(*(origin + 1.20 * y1), r"$\hat y_{br}=\hat y_1$", color=frame_red, fontsize=10.5, ha="right")
    ax.text(*(origin + 1.22 * zbr), r"$\hat z_{br}$", color=frame_green, fontsize=11, ha="right")
    ax.plot([origin[0], projection[0]], [origin[1], projection[1]], color="0.35", lw=1.3, ls=(0, (4, 3)))
    ax.plot([projection[0], point[0]], [projection[1], point[1]], color="0.45", lw=1.2, ls=(0, (3, 2)))
    ax.scatter(*projection, s=26, facecolor="white", edgecolor="black", zorder=6)
    ax.scatter(*point, s=34, color="black", zorder=6)
    ax.text(*(projection + np.array([0.12, -0.20])), r"$\mathbf{p}_{A/B}^{\mathrm{sag}}$", fontsize=10.5)
    ax.text(*(point + np.array([0.12, 0.03])), r"$\mathbf{p}_{A/B}$", fontsize=10.5)
    theta = np.arctan2(xbr[1], xbr[0])
    t = np.linspace(0, theta, 40)
    arc = origin + 0.82 * np.c_[np.cos(t), np.sin(t)]
    ax.plot(arc[:, 0], arc[:, 1], color=angle_color, lw=1.8)
    ax.text(*(origin + np.array([0.78, 0.28])), r"$-\theta_y$", color=angle_color, fontsize=11)
    ax.text(0.12, 3.48, "(B)", fontsize=12, weight="bold")
    ax.text(0.62, 3.48, r"Second rotation: $\hat x_{br}$ targets $\mathbf{p}_{A/B}$", fontsize=11)

    fig.subplots_adjust(left=0.03, right=0.98, bottom=0.02, top=0.99, hspace=0.08)
    fig.savefig(HERE / "xiFrameTransform.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(HERE / "xiFrameTransform.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    transform_pdf = HERE / "xiFrameTransform.pdf"
    reader = PdfReader(str(transform_pdf))
    writer = PdfWriter()
    writer.append_pages_from_reader(reader)
    alt_text = (
        "Two vertically stacked diagrams construct a bracket coordinate frame. "
        "Panel A rotates the intermediate x axis about the fixed space z axis until it points "
        "to the sagittal-plane projection of attachment point p A or p B. Panel B rotates "
        "about the current y axis until the final bracket x axis points to the full "
        "three-dimensional attachment point."
    )
    writer.add_metadata(
        {
            "/Title": "Two-rotation construction of the bracket frame",
            "/Subject": alt_text,
            "/AltText": alt_text,
        }
    )
    temp_pdf = HERE / "xiFrameTransform.metadata.pdf"
    with temp_pdf.open("wb") as stream:
        writer.write(stream)
    temp_pdf.replace(transform_pdf)


def _build_steele_composite_legacy():
    source_pdf = Path(
        r"C:\Users\Ben Bolen\Zotero\storage\4MQUGRSR\Steele et al. - 2022 - Experimental Verification of Kinematics and Kineti.pdf"
    )
    reader = PdfReader(str(source_pdf))
    source_images = reader.pages[3].images
    if len(source_images) < 2 or source_images[1].image.size != (843, 323):
        raise RuntimeError("Expected Steele et al. Figure 6 image was not found on PDF page 4")
    mechanism_image = source_images[1].image

    fig = plt.figure(figsize=(6.7, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=(1.18, 0.82), hspace=0.20)

    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 8.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(A)", loc="left", pad=6, fontsize=11.5, weight="bold")

    bone_fill = "#D9D9D9"
    bone_edge = "#555555"
    # Distal femur: shaft, flared metaphysis, and two condyles.
    ax.add_patch(
        patches.FancyBboxPatch(
            (4.35, 5.05),
            1.30,
            2.65,
            boxstyle="round,pad=0.08,rounding_size=0.28",
            facecolor=bone_fill,
            edgecolor=bone_edge,
            lw=1.4,
            zorder=0,
        )
    )
    ax.add_patch(
        patches.Polygon(
            [[4.30, 5.45], [3.55, 4.15], [4.05, 3.55], [5.95, 3.55], [6.45, 4.15], [5.70, 5.45]],
            closed=True,
            facecolor=bone_fill,
            edgecolor=bone_edge,
            lw=1.4,
            zorder=0,
        )
    )
    ax.add_patch(patches.Ellipse((4.20, 3.75), 1.25, 1.05, facecolor="#C8C8C8", edgecolor=bone_edge, lw=1.2, zorder=0))
    ax.add_patch(patches.Ellipse((5.80, 3.75), 1.25, 1.05, facecolor="#C8C8C8", edgecolor=bone_edge, lw=1.2, zorder=0))

    # Proximal tibia: plateau, metaphysis, and shaft.
    ax.add_patch(
        patches.Polygon(
            [[3.25, 2.45], [3.65, 3.05], [6.35, 3.05], [6.75, 2.45], [5.65, 1.78], [4.35, 1.78]],
            closed=True,
            facecolor=bone_fill,
            edgecolor=bone_edge,
            lw=1.4,
            zorder=0,
        )
    )
    ax.add_patch(
        patches.FancyBboxPatch(
            (4.42, 0.20),
            1.16,
            1.90,
            boxstyle="round,pad=0.06,rounding_size=0.20",
            facecolor=bone_fill,
            edgecolor=bone_edge,
            lw=1.4,
            zorder=0,
        )
    )

    f1 = np.array([3.82, 4.55])
    f2 = np.array([6.18, 4.40])
    t1 = np.array([3.90, 2.50])
    t2 = np.array([6.10, 2.42])
    # Link 1 and link 3 cross between femur and tibia; top/bottom links attach
    # the mechanism clearly to the two bone segments.
    ax.plot([f1[0], f2[0]], [f1[1], f2[1]], color="black", lw=7.0, solid_capstyle="round", zorder=3)
    ax.plot([t1[0], t2[0]], [t1[1], t2[1]], color="#555555", lw=7.0, solid_capstyle="round", zorder=3)
    ax.plot([f1[0], t2[0]], [f1[1], t2[1]], color="#CC6677", lw=5.5, solid_capstyle="round", zorder=2)
    ax.plot([f2[0], t1[0]], [f2[1], t1[1]], color=BLUE, lw=5.5, solid_capstyle="round", zorder=2)
    for point in (f1, f2, t1, t2):
        ax.add_patch(patches.Circle(point, 0.16, facecolor="white", edgecolor="black", lw=1.4, zorder=5))

    # Intersection of the crossed links marks the instantaneous center.
    A = np.column_stack((t2 - f1, -(t1 - f2)))
    alpha = np.linalg.solve(A, f2 - f1)[0]
    icr = f1 + alpha * (t2 - f1)
    ax.scatter(*icr, marker="+", s=180, linewidths=2.5, color="#AA3377", zorder=7)
    ax.text(*(icr + np.array([0.24, 0.06])), "ICR", color="#AA3377", fontsize=11, weight="bold")

    ax.annotate(
        "distal femur",
        xy=(5.10, 6.10),
        xytext=(7.25, 6.85),
        arrowprops=dict(arrowstyle="-", color="0.25", lw=1.1),
        fontsize=11,
        ha="left",
    )
    ax.annotate(
        "proximal tibia",
        xy=(5.05, 1.25),
        xytext=(7.20, 0.72),
        arrowprops=dict(arrowstyle="-", color="0.25", lw=1.1),
        fontsize=11,
        ha="left",
    )
    ax.annotate("femoral link", xy=(5.00, 4.48), xytext=(0.70, 5.10), arrowprops=dict(arrowstyle="-", lw=1.0), fontsize=10)
    ax.annotate("tibial link", xy=(5.00, 2.46), xytext=(0.82, 1.85), arrowprops=dict(arrowstyle="-", lw=1.0), fontsize=10)
    ax.text(7.32, 4.40, "crossed links", fontsize=10.5, color="0.25")

    ax2 = fig.add_subplot(gs[1])
    ax2.imshow(mechanism_image)
    ax2.axis("off")
    ax2.set_title("(B)", loc="left", pad=6, fontsize=11.5, weight="bold")

    fig.subplots_adjust(left=0.03, right=0.98, bottom=0.02, top=0.98)
    fig.savefig(HERE / "steeleknee.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(HERE / "steeleknee.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    alt_text = (
        "Panel A is a stylized sagittal schematic of the robotic knee with the distal femur and proximal tibia "
        "drawn as gray bone segments. Four pinned links form a crossed four-bar mechanism between them, and the "
        "intersection of the crossed links marks the instantaneous center of rotation. Panel B reproduces Steele "
        "et al. Figure 6, showing a sagittal cutaway of the knee linkage at full extension and an exploded view of "
        "the mechanical stop."
    )
    output_pdf = HERE / "steeleknee.pdf"
    output_reader = PdfReader(str(output_pdf))
    writer = PdfWriter()
    writer.clone_document_from_reader(output_reader)
    writer.add_metadata({"/Title": "Steele biomimetic knee mechanism", "/Subject": alt_text, "/AltText": alt_text})
    temp_pdf = HERE / "steeleknee.metadata.pdf"
    with temp_pdf.open("wb") as stream:
        writer.write(stream)
    temp_pdf.replace(output_pdf)


def build_steele_figure():
    """Reproduce only Steele et al. Figure 6, retaining its native (a)/(b) labels."""
    source_pdf = Path(
        r"C:\Users\Ben Bolen\Zotero\storage\4MQUGRSR\Steele et al. - 2022 - Experimental Verification of Kinematics and Kineti.pdf"
    )
    reader = PdfReader(str(source_pdf))
    source_images = reader.pages[3].images
    if len(source_images) < 2 or source_images[1].image.size != (843, 323):
        raise RuntimeError("Expected Steele et al. Figure 6 image was not found on PDF page 4")
    mechanism_image = source_images[1].image

    output_dir = HERE.parent / "Figures" / "Aim2"
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(7.2, 7.2 * 323 / 843))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.imshow(mechanism_image)
    ax.axis("off")
    fig.savefig(output_dir / "steeleknee.pdf", bbox_inches="tight", pad_inches=0.0)
    fig.savefig(output_dir / "steeleknee.png", dpi=300, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)

    alt_text = (
        "Two panels reproduced from Steele et al. Figure 6. Panel a is a sagittal cutaway "
        "of the biomimetic knee at full extension, showing the crossed four-bar linkage. "
        "Panel b is an exploded view of the knee's mechanical stop."
    )
    output_pdf = output_dir / "steeleknee.pdf"
    output_reader = PdfReader(str(output_pdf))
    writer = PdfWriter()
    writer.clone_document_from_reader(output_reader)
    writer.add_metadata(
        {
            "/Title": "Steele biomimetic knee mechanism",
            "/Subject": alt_text,
            "/AltText": alt_text,
        }
    )
    temp_pdf = output_dir / "steeleknee.metadata.pdf"
    with temp_pdf.open("wb") as stream:
        writer.write(stream)
    temp_pdf.replace(output_pdf)


if __name__ == "__main__":
    build_projection_figure()
    build_transform_figure()
    build_steele_figure()
    print("Wrote staged PDF/PNG figures to", HERE)
