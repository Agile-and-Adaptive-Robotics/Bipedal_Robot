# Expanded MN feedback and muscle maps

Open **Expanded_MN_Feedback_and_Muscle_Maps.pdf** for the two newly expanded examples. **RG_PF_MN_Full_Feedback.pdf** contains the full 16-page set, with these examples replacing the earlier PF-to-MN-only pages at pages 13 and 14. **Representative_MN_Layers_Expanded.pdf** contains the five MN examples alone.

The shared-PF sheet expands vastus medialis and vastus lateralis. Each muscle has its own Ia, II, Ib, excitatory and inhibitory interneurons, MN, and Renshaw cell. Their Renshaw cross-inhibition is retained. Four antagonist MN destinations are shown as a clearly marked bundle containing four distinct pools: semimembranosus, semitendinosus, and biceps femoris long and short heads. Every incoming line to that bundle represents four separate saved synapses. The four targets' own local loops are outside this crop. Bundling is permitted only after checking that the four connections on each line have identical conductance, reversal, threshold, and saturation parameters.

The biarticular sheet expands rectus femoris and semimembranosus. Each retains one muscle, one MN, one set of sensory maps, and one local reflex circuit. Hip and K&A PF projections converge on the same corresponding muscle MN. The actual saved Ia-IN and II-IN inhibition to the opposing MN, mutual Ia-IN inhibition, mutual Ib-IN inhibition, and Renshaw cross-inhibition are drawn.

Both sheets explicitly include these physical paths:

- MN → AD/map → muscle activation.
- Muscle-tendon length and velocity → spindle symbol → separate AD/maps → Ia and II inputs.
- Muscle force → tendon-force symbol → AD/map → Ib input.

AD is the circular adapter/map symbol, not an extra neuron. Brown arrows represent physical signals and current injection. Open triangles with the broad base on the target represent excitatory synapses; filled circles represent inhibitory synapses. Spindle and tendon symbols identify the encoded signals rather than adding new mechanistic receptor dynamics. The motor map uses the saved voltage-to-activation conversion, and the three sensory injection ports per expanded muscle are checked against the saved input-adapter entries. Source encoder and actuation policies are preserved in VERIFICATION.json.

Central afferent→RG/PF projections remain in the preceding RG/PF sheets. These two sheets expand local motor circuitry, descending PF inputs, and physical maps. The named w2lvar standing condition has zero effective PF-to-MN conductance; those source-defined routes remain gray and dashed. Other depicted neural connections are nonzero.

The source ledger is read as data only. Exact parity is checked for 82 depicted neural records; 40 directed physical signal paths are recorded separately. Every local and descending saved connection between the displayed neural cells is included. Geometry checks find no detached endpoints, lines hidden by nodes, or diagonal segments outside 45 degrees. Both rendered pages were visually inspected, including physical arrow contact with the custom muscle symbols and the preservation of dashed paths in PDFs.

Each sheet has editable SVG, PDF, PNG, neural-connection JSON, synapse CSV, physical-map CSV, and geometry-check files. No active model source or configuration was edited, and no circuit, training, optimization, or physics was run.
