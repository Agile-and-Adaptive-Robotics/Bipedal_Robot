import io

p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\50-database\components\sadb_figures_alt.txt"
add = """

--- App screenshots (text figures of record) ---

sadb_pivot_app.png - Screenshot of the SADb Explorer pivot view. A matrix heat
map with rows Ia, Ib, II, Cutaneous, Flexor reflex afferents, III/IV,
Mechanosensory, Nociceptive (and untagged) and columns Cat, Human, Rodent,
Vertebrates, Mammals, Invertebrates, and other animal tags; each cell shades
and counts the papers carrying both tags, dominated by the cat preparation.
App chrome around the matrix: corpus stats header, tab bar (Table, Pivot,
Bubble map, Search online), and the pivot dimension dropdowns.

sadb_clusters_app.png - Screenshot of the SADb Explorer bubble map view in
"Year x citations" layout, "Neurons" style, focused on Grillner 1985
("Neurobiological Bases of Rhythmic Motor Acts in Vertebrates", 1085
citations). Center of the map: the dark focus bubble. Green-triangle markers
connected by curved edges are the 57 in-corpus papers citing the focus; red
circles are the 13 papers the focus cites; grayed bubbles and edges are
second-degree connections, 555 papers within two degrees in all. Left panel
identifies the focused paper (title, author, year, citations, DOI
10.1126/science.3975635, animal tags Vertebrates/Mammals/Lamprey) with
References map / Cited-by map / Details actions; right panel lists the
connections ("Cited by 57 (excitatory input)" with entries such as Duysens et
al. 2004, Rossignol 2011, Taga 1991). X axis: publication year (ticks 1940 to
2020). Crumb at top reads "Focused: Grillner 1985 - ... 555 papers within 2
degree(s)".
"""
with io.open(p, "a", encoding="utf-8") as f:
    f.write(add)
print("alt text appended")
