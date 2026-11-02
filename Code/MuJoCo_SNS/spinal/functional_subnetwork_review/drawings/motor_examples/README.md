# Representative motor-neuron circuits

`RG_PF_and_Representative_MN_Layers.pdf` combines the existing 11 RG/PF pages with five new MN examples. The new examples begin on combined page 12. PDF bookmarks navigate both parts. `Representative_MN_Layers.pdf` contains just the five new examples.

Each example also has an editable SVG, a PNG preview, exact saved-connection JSON, synapse CSV, and geometry check beside it.

1. **Paired MN circuit:** the saved W2L108 hip circuit, with every saved connection among the displayed cells. This follows the supplied picture's paired arrangement. Native names are clarified in the caption: peripheral cells named “Ia IN” feed the central reciprocal cells named “Ia”. The saved RE-to-own-MN connection is excitatory; inhibitory symbols are retained on RE-to-Ia and reciprocal RE connections.
2. **Shared PF outputs:** three vastus MNs receive the same K&A extensor PF source, while two hamstring MNs receive the K&A flexor source. Each is a separate motor pool with a separate connection weight.
3. **Biarticular convergence:** rectus femoris receives Hip-F and K&A-E inputs; semimembranosus receives Hip-E and K&A-F inputs. Each muscle retains one MN pool. Both PF inputs are excitatory connections with separate conductances. They do not require duplicated hip and knee MNs or duplicated muscle reflex circuits.
4. **Multiple antagonist MNs:** one tibialis-anterior local Ia/II/Ib/Renshaw motif sends separate Ia-IN and II-IN inhibitory connections to three selected plantarflexor MNs. This expands representative local circuitry without repeating every muscle's motif.
5. **Active synergy convergence:** all three saved PF inputs to the selected rectus-femoris MN, alongside selected outputs of those channels to vastus medialis and soleus. Other PF-channel inputs to the latter two MNs are outside this view. These are NMF PF channels, not anatomical joint half-centers.

The w2lvar saved standing ledger has zero effective PF-to-MN conductance. Examples 2 and 3 therefore distinguish the builder's defined connection weights from the saved effective value. Gray dashed paths mean defined routing with effective g=0, not an active saved connection. The editable SVG and rendered PDFs both retain this distinction. No parameters were installed or changed to make the routes active.

Open triangles have their point on the incoming line and broad base on the target. Filled circles indicate inhibitory connections. Rounded wire bends and semicircular crossing bridges follow the supplied drawing preferences. Every depicted line terminates on a displayed object; each actual muscle MN remains distinct. Left-side source identities remain in the data files, while unnecessary side suffixes are omitted from figure labels.

`VERIFICATION.json` records source hashes and exact parity for 55 represented connection records. Geometry checks find no detached endpoints, hidden wires through nodes, or diagonal segments outside the permitted 45-degree direction. PNGs of all five pages were visually inspected. The PDF conversion explicitly preserves dashed lines because the SVG importer otherwise renders them solid.

No circuit, training, optimization, physics, native model execution, manuscript, or active source configuration was run or changed.
