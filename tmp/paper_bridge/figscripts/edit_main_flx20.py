import io
p = r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex'
s = io.open(p, encoding='utf-8', newline='').read()
orig = s

# 1. Preamble: subcaption after caption package
old1 = "\\usepackage{caption}\n"
new1 = "\\usepackage{caption}\n\\usepackage{subcaption}\n\\renewcommand\\thesubfigure{\\Alph{subfigure}}\n"
assert s.count(old1) == 1, 'anchor1'
s = s.replace(old1, new1)

# 2. Replace the stale flexor-pending sentence with a forward pointer
old2 = "The corresponding flexor measurements are also pending."
new2 = "The corresponding flexor validation appears below (Figure~\\ref{fig:Flx20mm})."
assert s.count(old2) == 1, 'anchor2'
s = s.replace(old2, new2)

# 3. Insert flexor paragraph + three-panel subfigure after the pressures paragraph
anchor3 = "The measured-versus-predicted torque comparison for the sweep follows once the reaction-point constants are entered and the adjoint values are recorded.\n"
assert s.count(anchor3) == 1, 'anchor3'
insert3 = anchor3 + """
Physical validation of the redesigned flexor route followed on the same test stand with the two \\qty{20}{\\mm} actuators of the design of record ($l_{0}=\\qty{56.2}{\\cm}$, Figure~\\ref{fig:Flx20mm}). Panel \\textbf{A} shows the target: the maximum human isometric flexion torque of the biceps femoris short head. Panel \\textbf{B} contrasts the two predictions of the same design at \\qty{620}{\\kPa} with the measurements: the original rigid-body prediction reaches a peak torque magnitude of \\qty{190.3}{\\newton\\meter}, while the corrected model, carrying the length offset, attachment compliance, and wrapping loss identified on the pinned knee, reaches \\qty{91.9}{\\newton\\meter}. Fifteen isometric points held from \\ang{-119} to \\ang{5} at hung loads of \\qtyrange{11.1}{106}{\\newton} and gauge pressures of \\qtyrange{226.5}{638}{\\kPa}, with the measured wrench carried to the moving knee center by the same adjoint transform as the identification data. Panel \\textbf{C} reads the verdict: the measured torque magnitude met or exceeded the human target at every point, with margins from \\qty{0.42}{\\newton\\meter} at \\ang{-105} to \\qty{17.0}{\\newton\\meter} at \\ang{-14}.

\\begin{figure}[H]
\\begin{center}
\\begin{subfigure}{0.32\\textwidth}
\\includegraphics[width=\\textwidth]{Figures/Flx20mm_humanTarget.eps}
\\caption{Maximum human isometric knee-flexion torque (biceps femoris short head).}\\label{fig:Flx20mm_A}
\\end{subfigure}\\hfill
\\begin{subfigure}{0.32\\textwidth}
\\includegraphics[width=\\textwidth]{Figures/Flx20mm_rigidVsCorrected.eps}
\\caption{Redesigned flexor route at \\qty{620}{\\kPa}: original rigid-body prediction, corrected prediction ($\\chi_{0}$--$\\chi_{3}$), and measured points.}\\label{fig:Flx20mm_B}
\\end{subfigure}\\hfill
\\begin{subfigure}{0.32\\textwidth}
\\includegraphics[width=\\textwidth]{Figures/Flx20mm_meetsHuman.eps}
\\caption{Corrected prediction and measured points against the human target.}\\label{fig:Flx20mm_C}
\\end{subfigure}
\\end{center}
\\caption{Physical validation of the redesigned flexor route (two $\\diameter$\\qty{20}{\\mm} BPAs, $l_{0}=\\qty{56.2}{\\cm}$) on the biomimetic knee. Knee angle is negative in flexion, and flexion torque is negative. Indigo filled circles are measured torque, with marker area scaling by gauge pressure; the black dotted line is the maximum human isometric torque; the magenta dashed line is the corrected prediction at \\qty{620}{\\kPa}; the gray dash-dot line in \\textbf{B} is the original rigid-body prediction at \\qty{620}{\\kPa}. The measured torque magnitude exceeds the human target at all fifteen points.}\\label{fig:Flx20mm}
\\end{figure}
"""
s = s.replace(anchor3, insert3)

# 4. Limitations: replace the stale sparse-validation clause
old4 = "the physical validation rests on three torque points near full extension, with the full sweep and the flexor measurements pending"
new4 = "the physical validation rests on sparse isometric holds, six points on the extensor cam route and fifteen on the flexor route, rather than dense sweeps of either configuration"
assert s.count(old4) == 1, 'anchor4'
s = s.replace(old4, new4)

assert s != orig
io.open(p, 'w', encoding='utf-8', newline='').write(s)
print('Main.tex edits applied: 4 anchors OK')
