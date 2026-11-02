"""Preserve final negative results and narrowed CAD subset qualifications."""
import json
from pathlib import Path
from update_family_cad_record import section
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')

def main():
 cad=ROOT/'cad_bpa';gait=ROOT/'gait_bpa';pair=json.loads((cad/'mount_subset/SLEEVE_PAIR_CLEARANCE.json').read_text())
 section(cad/'EXPERIMENT_RECORD.md','Simultaneous sleeve-interference check and adopted subset',
 'A subsequent mutual-sleeve check rejects the combined four-sleeve concept: the bilateral biceps-femoris-short envelopes overlap by 18.07 mm in the sampled historical motion even though each clears the CAD obstacles. The adopted simultaneous subset is therefore only the two vas_med sleeve envelopes. Their minimum mutual surface gap is 63.43 mm across the 17 sampled poses. The body-only four-placement report remains preserved as intermediate geometry; it must not be quoted as a simultaneous four-actuator mounting pass.\n\nThe final inspected figure mount_subset/CAD_SLEEVE_MOUNT_PROPOSAL.png/.pdf/.svg shows only this two-sleeve subset on actual CAD geometry. It is a placement proposal, not a pressure simulation or built robot. Fittings, Bowden cable route, brackets, attachment loads, mass, dynamics, and continuous swept motion remain unverified. SLEEVE_PAIR_CLEARANCE.json is the adopted subset ledger. This does not add a complete calibrated actuator assignment.\n\nThe first rendering attempt used a cleanup API absent in MuJoCo 2.3.7; the compatibility guard was corrected, and the final render was inspected. No numerical physics output was changed.')
 diag=json.loads((gait/'w2lvar_current_followup/CURRENT_BALANCE_DIAGNOSTIC.json').read_text())
 section(gait/'EXPERIMENT_RECORD.md','One-dimensional w2lvar current-balance diagnostic',
 'After the equal-budget comparison, one separate diagnostic holds standing bias at 0.8333333333333334 and aggregate Ib gain at zero, while setting shared PF→MN gain to 0.25. This tests the measured residual-inhibition explanation with no rewiring or broad sweep. Foot load rises to 0.63885 BW, but tilt reaches 13.080 degrees and sway 69.01 mm; the unchanged standing gate still fails. Mean upward rig force is 20.618% BW and downward force 3.021% BW. Thus a single global excitation multiplier trades load support against posture, rather than resolving the family mapping. No fresh pass, pressure transfer, or curriculum promotion is claimed. Exact case: w2lvar_current_followup/CURRENT_BALANCE_DIAGNOSTIC.json.\n\nThe prepared hybrid state audit, additional-family comparison, and CAD source/interface work are complete for this bounded follow-up. Sources and configurations are archived in followup_source_snapshot/MANIFEST.json; model/mesh and scientific-package dependencies remain external. The accepted 20 s syn6 guided-standing checkpoint is unchanged. No manuscript files or CHATGPT_REPORT.md were edited.')
 print('Final qualified follow-up record saved')

if __name__=='__main__':main()
