"""Lay out original OpenSim Desktop viewports, MuJoCo frames, and a Simulink diagram."""

import argparse
import io
from pathlib import Path
import subprocess
import tempfile

from PIL import Image
import pymupdf


def ensure_resolution(image):
    if image.height < 400:
        raise ValueError("The supplied viewport is too short for publication.")


def place_viewport(page, rect, source, crop_native_background=False):
    with Image.open(source) as original:
        if original.mode == "RGBA":
            bounds = original.getchannel("A").getbbox()
            if bounds is None:
                raise ValueError(f"The OpenSim snapshot is empty: {source}")
            image = original.crop(bounds)
        else:
            image = original.convert("RGB")
        if crop_native_background:
            # Keep the vertical field and remove only the far background at each side.
            image = image.crop((round(.245 * image.width), 0, round(.755 * image.width), image.height))
        ensure_resolution(image)
        contents = io.BytesIO()
        image.save(contents, "PNG", optimize=True)
        page.insert_image(rect, stream=contents.getvalue(), keep_proportion=True)


def native_diagram_rect(page):
    bounds = pymupdf.Rect()
    for _, coordinates in page.get_bboxlog():
        rect = pymupdf.Rect(coordinates)
        if rect.get_area() < .75 * page.rect.get_area():
            bounds |= rect
    # Simulink prints the diagram onto a white landscape page.
    bounds += (-3, -3, 3, 3)
    return bounds & page.rect


def panel_letters(page, positions, font):
    for i, x in enumerate(positions):
        page.insert_text((x, 13), f"({chr(65+i)})", fontfile=str(font),
                         fontname="PanelArialBold", fontsize=13, color=(0,0,0))


def save(doc, path):
    doc.subset_fonts()
    doc.save(path, garbage=4, deflate=True)
    doc.close()


def produce_native(output, components, model_root, python, font, azimuth, elevation):
    script = components / "native_mujoco.py"
    stem = "gait2392_"
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        for kind, name in [("human", "simbody"), ("robot", "robot")]:
            base = model_root / "mjc" / (stem + name)
            cmd = [str(python), "-B", str(script), str(base / (stem + name + "_cvt3.xml")),
                   str(work / (kind + ".png")), "--azimuth", str(azimuth), "--elevation", str(elevation)]
            if kind == "robot":
                cmd.extend(["--bridge", str(components / "robot_bridge_readback.json.gz")])
            subprocess.run(cmd, check=True)

        human = pymupdf.open()
        page = human.new_page(width=191, height=150)
        panel_letters(page, (5, 102), font)
        place_viewport(page, pymupdf.Rect(3, 17, 92, 149), components / "opensim_human.png")
        place_viewport(page, pymupdf.Rect(100, 17, 189, 149), work / "human.png", True)
        save(human, output / (stem + "human_models.pdf"))

        robot = pymupdf.open()
        page = robot.new_page(width=390, height=150)
        panel_letters(page, (4, 85), font)
        place_viewport(page, pymupdf.Rect(2, 17, 76, 149), components / "opensim_robot.png")
        place_viewport(page, pymupdf.Rect(83, 17, 143, 149), work / "robot.png", True)
        with pymupdf.open(components / "robot_interface.pdf") as diagram:
            page.show_pdf_page(pymupdf.Rect(148, 17, 387, 149), diagram, 0,
                               clip=native_diagram_rect(diagram[0]), keep_proportion=True)
        save(robot, output / (stem + "robot_models.pdf"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--mujoco-python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--font", type=Path, default=Path("C:/Windows/Fonts/arialbd.ttf"))
    parser.add_argument("--components", type=Path, default=Path(__file__).parent)
    parser.add_argument("--azimuth", type=float, default=135)
    parser.add_argument("--elevation", type=float, default=-12)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    produce_native(args.output, args.components,
                   args.repo / "Solid_Models/OpenSim/Gait2392_Robotbody",
                   args.mujoco_python, args.font, args.azimuth, args.elevation)
