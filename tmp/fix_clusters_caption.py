# Swap the fig:sadb_clusters lead-in + caption to the focused Grillner-1985 app view.
# The figure file was replaced with the focused bubble-map screenshot on 2026-10-19.
import io

p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\50-database.tex"
t = io.open(p, encoding="utf-8").read()

old_lead = (
    "Figure~\\ref{fig:sadb_clusters} summarizes the same corpus at the topic-cluster level: "
    "the seven clusters the curated corpus retains, and the satellite clusters the prune removed in full."
)
new_lead = (
    "Figure~\\ref{fig:sadb_clusters} shows the same corpus as the packaged citation graph: "
    "the bubble map draws every paper as a neuron, its soma sized by citation count, and "
    "focusing a paper colors its direct citation neighborhood by direction."
)

old_cap = (
    "\\caption[The SADb Explorer bubble map view]{The SADb Explorer's bubble map view, the "
    "packaged citation graph of the full audited corpus: every paper is drawn as a neuron, its soma "
    "sized by citation count, with citation links rendered as excitatory and inhibitory synapses and "
    "cluster labels grouping the topic communities. Clicking a bubble focuses its citation neighborhood "
    "out to a chosen degree, with direct connections colored by direction. The view shows the curated "
    "core clusters, from the afferent-control classics and spinal-circuit modeling to neuromechanics and "
    "robot walking, alongside the satellite clusters the prune later removed.}"
)
new_cap = (
    "\\caption[The SADb Explorer bubble map view, focused on one paper]{The SADb Explorer's bubble "
    "map view, focused on a classic of the corpus (Grillner, 1985). Every paper is drawn as a neuron, "
    "its soma sized by citation count; in the focused view, direct connections are colored by direction "
    "(triangles mark the 57 in-corpus papers citing the focus, circles the 13 it cites) and second-degree "
    "papers are grayed, here 555 papers within two degrees. The left panel identifies the focused paper "
    "and the right panel lists its connections; the axes give publication year and citation count.}"
)

assert old_lead in t, "lead-in not found"
assert old_cap in t, "caption not found"
t = t.replace(old_lead, new_lead).replace(old_cap, new_cap)
io.open(p, "w", encoding="utf-8", newline="\n").write(t)
print("caption + lead-in updated")
