"""Apply the manuscript's panel callout style and the author's caption edits."""
from pathlib import Path
import re


def end_of_group(text, start, left, right):
    depth = 1
    for index in range(start + 1, len(text)):
        if text[index - 1] == "\\":
            continue
        if text[index] == left:
            depth += 1
        elif text[index] == right:
            depth -= 1
            if depth == 0:
                return index
    raise ValueError(f"Unclosed LaTeX group at {start}")


def figure_captions(text):
    floats = re.finditer(r"\\begin\{((?:sideways)?figure\*?)\}([\s\S]*?)\\end\{\1\}", text)
    for figure in floats:
        for command in re.finditer(r"\\caption\b\s*(?=[\[{])", figure[2]):
            start = figure.start(2) + command.end()
            if text[start] == "[":
                start = end_of_group(text, start, "[", "]") + 1
                while text[start].isspace():
                    start += 1
            assert text[start] == "{"
            stop = end_of_group(text, start, "{", "}")
            yield start + 1, stop, text[start+1:stop]


edits = {
    "10-introduction": [
        ("Dissertation roadmap. The solid arrows show relationships among the chapters. Dashed connections indicate future integration in the physical robot and the intended return of biological insight.",
         "Relationships among the dissertation studies and their proposed integration in the physical robot.")],
    "30-torque": [
        (" The legend distinguishes measurements, the Hunt model, the improved BPA model, and the hybrid calculation.", ""),
        (" Circles show measurements, dashed curves the improved BPA model, and solid magenta curves the optimized model.", ""),
        (" The legend applies to all panels.", ""),
        (r" The legend in \textbf{(B)} identifies the measurements, predictions, and human reference.", "")],
    "60-sns": [
        ("seven suspended syn6 cycles (light blue) and their mean (solid blue)", "seven suspended syn6 cycles and their mean"),
        ("OpenSim Subject01 example (dashed)", "OpenSim Subject01 example"),
        (r"\citet{ong_predicting_2019} (dotted)", r"\citet{ong_predicting_2019}"),
        ("Dark bars below each leg mark stance in each reference, and gray bars mark swing.", "The bars below each leg mark stance and swing in the walking references."),
        (" Solid lines show the left leg and dashed lines the right.", "")],
    "93-AppendixB": [
        ("(row A1)", r"(row \textbf{A1})"), ("(row B1)", r"(row \textbf{B1})"), ("(row C1)", r"(row \textbf{C1})")],
    "94-AppendixC": [
        ("Black dashed curves show the maximum-isometric human reference.", "The human reference gives maximum isometric torque."),
        ("poses A--F", r"poses \textbf{A--F}"), ("poses G and H", r"poses \textbf{G} and \textbf{H}"),
        ("Blue lines show the route; colored outlines show the clearance envelopes. Gray silhouettes show the registered OpenSim anatomical reference; no foot is drawn.",
         "The route and clearance envelopes are registered to the OpenSim femur and tibia.")],
    "95-AppendixD": [
        ("Solid blue curves show actual transmission; dashed black curves show the intended primitive, and dotted lines mark the 0.2 activity threshold.",
         "The dotted lines mark the 0.2 activity threshold."),
        ("Blue curves and bands show the model mean and one standard deviation across ten cycles per leg; black dashed curves show the reference.",
         "The model curves summarize ten cycles per leg. Shading spans one standard deviation about the mean."),
        (" Right is solid blue; left is dashed orange.", "")],
}


def revise_captions(chapter, text):
    for a, b, caption in reversed(list(figure_captions(text))):
        for old, new in edits.get(chapter, ()):
            caption = caption.replace(old, new)
        # The negative lookbehind preserves pre-existing bold labels. Parentheses
        # around an acronym or an equation are not panel references.
        caption = re.sub(r"(?<!\\textbf\{)(\([A-H](?:[1-9])?(?:(?:, |--)[A-H][1-9]?)*\))",
                         r"\\textbf{\1}", caption)
        text = text[:a] + caption + text[b:]
    return text


if __name__ == "__main__":
    chapters = Path(__file__).resolve().parents[2] / "Final/chapters"
    for source in chapters.glob("*.tex"):
        if source.stem == "20-muscles":
            continue  # Published chapter reviewed separately against its original.
        raw = source.read_bytes()
        new = revise_captions(source.stem, raw.decode("utf8")).encode("utf8")
        if new != raw:
            staged = source.with_suffix(".caption-edit")
            staged.write_bytes(new)
            staged.replace(source)
