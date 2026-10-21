import io

F = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
t = io.open(F, encoding="utf-8").read()

a1 = (r"The direct extensor route exceeds the human envelope in magnitude at every "
      r"angle, with zero shortfall points: across the design range the predicted "
      r"two-BPA extensor torque spans \qty{35.45}{\newton\cdot\meter} to "
      r"\qty{109.4}{\newton\cdot\meter} against the human-scaled requirement of "
      r"\qty{20.25}{\newton\cdot\meter} to \qty{50.32}{\newton\cdot\meter}, and the "
      r"required 5\,\% overshoot is consumed as the thinnest margin exactly as "
      r"designed, at a minimum margin of $+5.02$\,\% at the worst-case angle of "
      r"$-20.2^{\circ}$.")

fig1 = r"""

\begin{figure}[p]
\centering
\includegraphics[width=0.98\textwidth]{ExtensorRoute_redesigned.pdf}
\caption[Route geometry of the redesigned extensor route]{Route geometry of the direct extensor route redesign, drawn in the femur frame. Eight pose panels sweep the design range from full flexion ($\theta_\text{k} = \ang{-120}$) through the frame immediately before each of the six contact-elimination events to full extension ($\theta_\text{k} = \ang{10}$); knee angle is negative in flexion. Each panel draws the active route over the femur-cylinder, femur-wall, condyle, and tibia clearance envelopes of the route-elimination rule, with the local bend-radius checks dashed. The square marker is the femur-side origin and the diamond marker is the tibial insertion.}
\label{fig:ExtensorRoute_redesigned}
\end{figure}
"""

a2 = (r"On the knee-flexor route, the pulley variant converged to a collision-feasible "
      r"layout with zero bound violations and a peak predicted flexion torque of "
      r"\qty{85.56}{\newton\cdot\meter} on the insertion-side measure, the torque "
      r"applied at the insertion after the tackle's force loss; "
      r"Section~\ref{sec:disc_pulley} quantifies what the redirect adds on the "
      r"knee-axis measure and the off-axis torque it bounds.")

fig2 = r"""

\begin{figure}[p]
\centering
\includegraphics[width=0.98\textwidth]{PulleyTorqueSurfaces_20261020-eps-converted-to.pdf}\par\medskip
\includegraphics[width=0.98\textwidth]{PulleyMarginSurfaces_20261020-eps-converted-to.pdf}
\caption[Pulley-position sweep surfaces for the reeved flexor route]{Pulley-position sweep for the reeved knee-flexor route. Top: the predicted torque surface over the two pulley positions along the tendon run; bottom: the corresponding torque-margin surface against the human-scaled flexor requirement, on the same position grid. The sweeps show where the reeving freedom pays, the feasible ridges of pulley placement, and where the margin collapses as the geometry approaches the clearance constraints.}
\label{fig:PulleySurfaces}
\end{figure}
"""

if a1 in t:
    t = t.replace(a1, a1 + fig1)
    print("ext figure inserted")
else:
    print("EXT ANCHOR NOT FOUND")

if a2 in t:
    t = t.replace(a2, a2 + fig2)
    print("pulley figures inserted")
else:
    print("PULLEY ANCHOR NOT FOUND")

io.open(F, "w", encoding="utf-8").write(t)
