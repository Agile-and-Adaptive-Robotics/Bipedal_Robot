from pathlib import Path
import re
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
p=C/'chapters/60-sns.tex';s=p.read_text()
for old,new in [
('These were integration requirements of this project, not capabilities absent from computational neuroscience as a whole.','These requirements concerned integration for this robot project.'),
('The relevant distinction is the added actuator and mechanical interfaces, rather than graphical neural modeling itself.','The added components connect neural signals to the actuator and mechanical interfaces.'),
('The representative results below establish those capabilities.',''),
('This comparison demonstrates the response of the implemented feedback pathway. It does not validate the provisional mechanical parameters against an experiment or establish quantitative control of the physical knee.','The comparison checks the implemented feedback response; experimental calibration of the mechanical parameters remains necessary.'),
('These results verify agreement for this circuit and protocol, rather than interchangeability of arbitrary controllers or simulators.','The agreement is specific to this circuit and input protocol.'),
("The agreement of that negative result matters: this network's alternating activity depends on the specified sensory drive. It should not be interpreted as an autonomous oscillator merely because it alternates under the full protocol.","The alternating activity therefore depends on the specified sensory drive."),
('These tests establish that commands and measurements pass through the interface; they do not establish successful gait control. In particular,','Commands and measurements passed through the interface. However,'),
('The supported robot evaluations produced stepping-like motion and, in some configurations, forward progression. They did not demonstrate harness-free walking. Removing a tonic neural drive in a supported model does not remove the mechanical support. Spiking-neuron comparisons likewise preserved a continuous motor-activation interface but changed rhythm timing and, in some runs, limb coordination. The detailed comparisons are reported in Appendix~\\ref{app:simulation_verification}; they do not support a general claim that changing neural coding leaves movement unchanged.',r'The supported robot evaluations produced stepping-like motion and, in some configurations, forward progression, but no harness-free walking. Reducing tonic neural drive did not remove mechanical support. Spiking-neuron comparisons preserved a continuous motor-activation interface while changing rhythm timing and, in some runs, limb coordination (Appendix~\ref{app:simulation_verification}).'),
('They are not prerequisites for using the demonstrated reflex, conversion, and cross-implementation tools.','')]:
    s=s.replace(old,new)
p.write_text(s,encoding='utf-8')
def modify_table(s,label,old,new):
    for m in re.finditer(r'\\begin\{table\}.*?\\end\{table\}',s,re.S):
        if r'\label{'+label+'}' in m.group():
            block=m.group();assert old in block,(label,old)
            block=block.replace(old,new).replace(r'\end{tabular}',r'\end{tabularx}')
            return s[:m.start()]+block+s[m.end():]
    raise ValueError(label)
p=C/'chapters/95-AppendixD.tex';s=p.read_text()
s=modify_table(s,'tab:walker_ablations',r'\begin{tabular}{@{}llr@{}}',r'\small\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.20\textwidth}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}p{0.24\textwidth}@{}}')
s=modify_table(s,'tab:spike_knee',r'\begin{tabular}{@{}lcc@{}}',r'\small\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.42\textwidth}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X@{}}')
s=s.replace('Untuned spiking mirror versus its non-spiking twin at the same winner parameters, ground evaluation (16 s). The cycle knee minimum and stance fractions are computed over detected gait cycles; the kinematic score is the cycle-shape error metric of the production pipeline (0 denotes perfect agreement; higher scores indicate closer agreement). The spiking mirror carries mapped parameters only, with zero re-tuning.','Supported ground evaluation at the same selected parameters (16 s), with no spiking re-tuning. Knee minima and contact duty are computed over detected cycles. The kinematic score is a cycle-shape error metric: zero denotes perfect agreement, and less negative scores indicate closer agreement.')
s=s.replace('committed baseline versus spiking twin','graded baseline versus spiking twin')
p.write_text(s,encoding='utf-8')
p=C/'chapters/96-AppendixE.tex';s=p.read_text().replace('p{0.08\\textwidth}','p{0.14\\textwidth}').replace('p{0.22\\textwidth}','p{0.16\\textwidth}');p.write_text(s,encoding='utf-8')
p=C/'chapters/94-AppendixC.tex';s=p.read_text()
s=modify_table(s,'tab:bipulley_batch',r'\begin{tabular}{@{}lll@{}}',r'\small\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.22\textwidth}>{\raggedright\arraybackslash}p{0.32\textwidth}>{\raggedright\arraybackslash}X@{}}')
p.write_text(s,encoding='utf-8')
print('Adjusted new appendix tables and tightened repeated qualifications')
