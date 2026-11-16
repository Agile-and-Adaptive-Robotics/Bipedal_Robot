"""Apply independently verified historical scope and chronology corrections."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
CH=HERE.parent.parent/'Final/chapters'
fixes={
'60-sns.tex':[
('bridge verification also tests a different stated rate combination.',
 r'the bridge two-rate test uses a \qty{1}{\milli\second} command-source interval and a \qty{5}{\milli\second} mechanical/sensor interval.'),
('The five tested sensor signals were bit-exact over the scored 80-step window.',
 r'The five tested sensor arrays were bit-exact over 80 steps from \qty{0.005}{\second} to \qty{0.400}{\second}, excluding the initial sample.'),
(r'For the contact-driven W2L rhythm core, identical neuron definitions, synapses, parameters, drive, and contact inputs were supplied to NumPy and Simulink implementations. Burst periods matched at \qty{1.000}{\second}, all ten expected bursts were reproduced, and voltage correlation was at least 0.99996 (Figure~\ref{fig:w2l_crossengine}). Removing the phasic heel input caused both realizations to latch. This verifies the tested circuit calculation and its dependence on the specified inputs, before introducing a different mechanical feedback history.',
 r'The contact-driven W2L comparison tested a representative 26-neuron Simulink rhythm core against the corresponding signals of the 95-neuron NumPy network. The core used corresponding neuron parameters and shared pathways under identical drive and prescribed heel inputs. Over the post-startup \qtyrange{2}{12}{\second} window, burst periods matched at \qty{1.000}{\second} and all ten expected bursts were reproduced. Best-lag correlation of the four shared RG voltage traces was at least 0.99996, allowing a \qty{2}{\milli\second} relative alignment lag (Figure~\ref{fig:w2l_crossengine}). Removing the phasic heel input caused both realizations to latch. This verifies the tested shared rhythm core and input dependence, rather than complete-network or unshifted sample agreement.'),
('The result identifies a contribution to contact coordination, rather than a requirement for all rhythmic movement (Figure~',
 r'The retained evaluator also specifies a common \qty{40}{\newton} forward pulse from \qty{8}{\second} to \qty{8.15}{\second} in both conditions; an unforced baseline was not established by this pair. The result identifies a contribution to contact coordination under those matched inputs, rather than a requirement for all rhythmic movement (Figure~'),
(r'None of the program\'s recorded results establishes harness-free walking.',
 r'The coupled controller-transfer demonstrations reported here do not establish harness-free walking.'),
],
'95-AppendixD.tex':[
('MuJoCo vs.\\ OpenSim muscle force along the subject01 walk at\n  matched activations.',
 r'MuJoCo vs.\ OpenSim right-side muscle force at matched activations along the subject01 reference interval, \qtyrange{0.5}{2.5}{\second} (121 poses).'),
(r'The contact-driven W2L rhythm-core comparison uses identical neuron parameters, connections, descending drive, and contact inputs in NumPy and Simulink. Burst periods matched at \qty{1.000}{\second}, all ten expected bursts were reproduced, and voltage-trace correlation was at least 0.99996 (Figure~\ref{fig:w2l_crossengine}). Removing the phasic heel inputs caused latching in both implementations. This result verifies the specified input-dependent rhythm core; it does not demonstrate that the circuit oscillates autonomously or that every larger network agrees across implementations.',
 r'The contact-driven W2L comparison maps corresponding neuron parameters and shared rhythm pathways from the 95-neuron NumPy network into a representative 26-neuron Simulink core. Drive and prescribed heel inputs match. Over the post-startup \qtyrange{2}{12}{\second} interval, burst periods matched at \qty{1.000}{\second} and all ten expected bursts were reproduced. The correlation comparison uses 5,000 matching \qty{2}{\milli\second} samples through \qty{11.998}{\second}; best-lag correlations of the four RG traces were at least 0.99996 at a \qty{2}{\milli\second} relative lag (Figure~\ref{fig:w2l_crossengine}). Removing the phasic heel inputs caused latching in both implementations. The shared-signal agreement does not establish identical complete networks, autonomous oscillation, or sample-exact transfer.'),
('NumPy and Simulink implementations of the contact-driven W2L rhythm core under identical drive and contact inputs. Agreement concerns the displayed neural circuit and prescribed inputs, rather than independent closed-loop mechanical trajectories.',
 r'Corresponding W2L rhythm signals in the full NumPy network and representative Simulink core under identical drive and prescribed heel inputs. The reported correlation is a best-lag measure with \qty{2}{\milli\second} relative alignment; it is not unshifted sample agreement or a comparison of independent closed-loop mechanical trajectories.'),
(r'The saved human-model controller reproduced its cycle-shape score of $-145.07$ in an independent process.',
 r'The saved human-model controller reproduced its composite kinematic score of $-145.07$ in an independent process. The analysis uses the post-startup \qtyrange{2}{16}{\second} recording; the score combines several kinematic terms and is not solely a cycle-shape error.'),
(r', and fitted speed of \qty{0.0918}{\meter\per\second}.',
 r'. The fitted speed of \qty{0.0918}{\meter\per\second} uses the additional \qtyrange{8}{30}{\second} fitting window.'),
(r'its mean vertical force still supplied approximately \qty{16.76}{\percent} of body weight.',
 r'its mean vertical force still supplied approximately \qty{16.75}{\percent}, using the legacy \qty{411}{\newton} weight denominator over \qtyrange{4.5}{30}{\second}.'),
('Mild and strong drive reductions, forward and backward \\qty{40}{\\newton} pushes, and separate Ia, Ib, Renshaw, crossed-heel, and heel-cascade deletions retained the supported-stepping criterion.',
 r'The control configuration inherits a \qty{40}{\newton} forward pulse from \qty{8}{\second} to \qty{8.15}{\second} in the retained evaluator. Its explicitly named forward-push recording is identical to the control in every saved numerical array, so it is not an independent push-versus-unforced comparison. Mild and strong drive reductions, reversal to a backward \qty{40}{\newton} pulse, and separate Ia, Ib, Renshaw, crossed-heel, and heel-cascade deletions retained the supported-stepping criterion.'),
('The Di Russo source mechanical model initialized successfully, but its native custom controller runtime was unavailable for a baseline replay.',
 r'The custom SCONE runtime for the model of \citet{dirusso_investigating_2023} was subsequently built, and a \qty{10}{\second} replay in the original mechanical model reproduced the native baseline, including an exact saved fresh-process repeat. This verifies execution of that reference stack; it does not validate transfer to Gait2392 or establish physiological validity.'),
('The retained source inventory and rejected-run records identify requirements for a future replication; the unsuccessful search is not evidence of walking.',
 'The native baseline and rejected adaptation are distinct results: reproducing the original stack supplies a reference for transfer, while the unsuccessful search remains rejected as walking.'),
('during the post-startup evaluation window, whether or not cycles were detected.',
 r'over $2\leq t<16$\,s, whether or not cycles were detected.'),
('sustained unilateral firing; one side at', 'strongly asymmetric firing; one side at'),
]}
# Match an ordinary apostrophe in the manuscript.
fixes['60-sns.tex'][-1]=( "None of the program's recorded results establishes harness-free walking.",
 'The coupled controller-transfer demonstrations reported here do not establish harness-free walking.')
for name,pairs in fixes.items():
    p=CH/name;text=p.read_text(encoding='utf-8')
    for old,new in pairs:
        if old not in text:
            if new in text:continue
            raise ValueError((name,old))
        if text.count(old)!=1:raise ValueError(('ambiguous',name,old))
        text=text.replace(old,new,1)
    p.write_text(text,encoding='utf-8',newline='')
print('Applied historical evidence findings to the two canonical files.')
