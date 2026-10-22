from pathlib import Path
import json,hashlib,shutil,difflib,ast
H=Path(__file__).resolve().parent;R=H.parent.parent
D=R/'Documentation/Reports and Papers/Dissertation';F=D/'Figures';P=D/'ProofFinal'
B=Path((H/'build_path.txt').read_text());V=json.loads((H/'validation.json').read_text())
OUT=P/'build_schematic_20261022';OUT.mkdir(exist_ok=True)
for src in B.iterdir():
    if src.is_file():shutil.copy2(src,OUT/src.name)
if (B/'chapters').exists():shutil.copytree(B/'chapters',OUT/'chapters',dirs_exist_ok=True)
N=D/'Notes/figure_polish_20261022';N.mkdir(exist_ok=True)
shutil.copy2(H/'validation.json',N/'validation.json')
shutil.copy2(H/'protected_hashes.json',N/'protected_hashes.json')
baseline=json.loads((H/'baseline_hashes.json').read_text());changes=[]
for path,old in baseline.items():
    p=Path(path)
    if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()!=old:changes.append(str(p.relative_to(R)))
(N/'changed_existing_files.json').write_text(json.dumps(changes,indent=2))
diff=[]
for n in ['10-introduction','60-sns','70-overall']:
    p=P/'chapters'/f'{n}.tex';q=H/'before'/p.relative_to(R)
    if q.exists():diff.extend(difflib.unified_diff(q.read_text().splitlines(True),p.read_text().splitlines(True),fromfile='before/'+n,tofile='after/'+n))
(N/'source_changes.diff').write_text(''.join(diff),encoding='utf-8')
for n in ['make_pptx_figures.py','export_pptx.py']:ast.parse((F/n).read_text(encoding='utf-8'))
protected=json.loads((H/'protected_hashes.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
assert V['pages']==298 and not any(V[k] for k in ['undefined_refs','undefined_citations','missing_files','protected_changed'])
for name in ['pipeline','bpa_construction','bpa_exploded_cutaway']:
    shutil.copy2(H/'renders'/f'{name}.png',N/f'{name}_preview.png')
report=f'''# Schematic and flowchart polish

## Scope and directives

Read CHATGPT_REPORT.md, Notes/new_session_prompt_20261022.md, the dissertation skill, and the required campaign context. Improved 11 editable components across nine numbered figures. No commits, pushes, or Overleaf uploads were made.

Chapter 3 (isometric force, source 20-muscles.tex), including every MATLAB figure, was protected. All {len(protected)} protected file hashes match the pre-task snapshot. Experimental data, reported results, bibliography, controller code, and simulations were not changed. Other pre-existing working-tree changes were preserved.

## Figures

| Figure | Components | Main change |
|---|---|---|
| 1.1 | roadmap | Consistent hierarchy, alignment, and distinction between developed work and future integration. |
| 2.2 | actuator_technology_comparison | Clearer mechanism drawings, consistent typography, and outward BPA pressure arrows. |
| 6.1 | pipeline | Readable simulation loop, hardware path, local sensors, and future return connections. |
| 6.6 | sns_library_overview; simulink_neuron_detail | Consistent library blocks and readable functional neuron schematic, verified against the saved Simulink implementation. |
| 6.7 | KneeReflex_circuit | Clear excitatory/inhibitory terminals, sensory grouping, descending drive, mechanical feedback, and separated labels. |
| 6.13 | sns_curriculum | Aligned stage rows, explicit objectives, and retained harness conditions. |
| 7.1 | design_workflow | Clear decision routes and verification before committing hardware. |
| 7.2 | bpa_construction; bpa_exploded_cutaway | Six-layer FEA schematic, contact interfaces, pressure and boundary conditions; corrected cutaway symmetry-plane callout and editable labels. |
| 7.3 | sensorium_hierarchy | Consistent control hierarchy, explicit proposed modules, and afferent returns to rhythm and motor layers. |

Native editable PowerPoint files, EPS exports, and alt text are under each chapter's Figures/components folder. Chapter-level PDFs are used by LaTeX. The cutaway rendering is an illustrative raster with native editable PowerPoint annotations; its prompt is preserved beside the source. The original Simulink detail export is preserved as 60-sns/components/simulink_neuron_detail_native_export.pdf.

Figures/make_pptx_figures.py and Figures/export_pptx.py reproduce the native sources and exports using the saved manifest. Other plots, photographs, borrowed circuit illustrations, and the protected Chapter 3 figures were retained. Existing unused graphics were preserved rather than silently deleted.

## Source edits and compliance

Only three chapter source files were adjusted for these schematic changes. The exact pre-task-to-final differences are in source_changes.diff.

- 10-introduction.tex:50: roadmap caption scopes solid blocks to work developed or used and dashed blocks to future integration.
- 60-sns.tex:118–119: neuron inclusion uses readable full width; caption describes the functional schematic, current summation, integration, inverse-slope mapping, saturation, and voltage feedback. The nearby introduction to the panel is updated consistently.
- 60-sns.tex:413: curriculum caption states suspended air stepping and retained ground support, replacing the misleading shading explanation.
- 70-overall.tex:40 and the following caption: uses the annotated cutaway PDF and distinguishes longitudinal symmetry from a movable, optionally loaded axial boundary. It explicitly discloses that the middle rubber layer is separated in panel A only.

No silent deletion of numerical results or alteration of scientific data was made. Redundant caption wording was shortened to keep both BPA panels and their caption together; the technical qualifications remain explicit.

## Verification

- Full MiKTeX build: pdflatex, bibtex with absolute bibliography paths, then two pdflatex passes, all exit 0.
- Final PDF: 298 pages at ProofFinal/build_schematic_20261022/main.pdf.
- Undefined references: 0; undefined citations: 0; missing files: 0.
- All nine changed figure pages rendered and inspected; no out-of-page text or clipped captions.
- Effective editable figure text is at least 10 pt (minimum measured 10.368 pt; other components 10.8 pt or larger).
- Independent technical review: PASS for all 11 components. The neuron equations and topology were checked against the saved library, including inverse-slope gain and saturation. The BPA callouts identify retained material and the correct symmetry plane.
- Independent visual review: all 11 components pass; final rebuilt reflex-page verification is recorded separately below.
- Meaning survives grayscale through explicit labels, placement, direction, and line styles. Crossings without connections have clear gaps where needed.
- All {len(protected)} protected Chapter 3 files remain byte-identical.

The app's built-in compiler was attempted, but its Tectonic backend rejects pre-existing EPS inputs in protected Chapter 3. They were not converted or modified. The canonical MiKTeX PDF succeeds. Existing unrelated table/appendix overflow warnings remain outside this schematic-only scope; this is not a claim that the entire dissertation has no layout warnings.

## Review and commit scope

Review the final PDF and the three preview PNGs in this folder. Commit only the intended schematic PDFs, native sources/exports/alt text, generator and manifest, three chapter adjustments, and this report folder. Keep unrelated pre-existing dissertation, plot, hardware, and CAD changes separate. The baseline copies and detailed working renders are retained in tmp/schematic_polish_20261022/before and its sibling folders. No existing top-level dissertation PDF or build_r8 was replaced.
'''
(N/'REPORT.md').write_text(report,encoding='utf-8')
print('PACKAGED',OUT,'REPORT',N,'PROTECTED',len(protected),'CHANGED_BASELINE',len(changes))
