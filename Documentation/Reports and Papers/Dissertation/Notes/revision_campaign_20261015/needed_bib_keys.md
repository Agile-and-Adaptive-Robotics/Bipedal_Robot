# needed_bib_keys.md — bibliography keys requested by chapter writers, 2026-10-15 campaign

Writers never edit thesis.bib. Append `key + full reference` here; the steward merges
them into thesis.bib. Verified against the AARL Zotero group library (735051) where
noted.

## Requested by 40-ankle.tex (Writer — Ankle chapter, 2026-10-15)

- `peterka_2002` — cited in 40-ankle.tex as \citet (protocol source for the
  platform-oscillation balance tests). In the AARL Zotero group (item UIUWTJ4K;
  duplicate record PW7XTG4G).

```bibtex
@article{peterka_2002,
	title = {Sensorimotor integration in human postural control},
	volume = {88},
	number = {3},
	pages = {1097--1118},
	journal = {Journal of Neurophysiology},
	author = {Peterka, R. J.},
	year = {2002},
	doi = {10.1152/jn.2002.88.3.1097},
}
```

## Requested by 15-background.tex (Writer — Background chapter, 2026-10-15)

New keys cited in the revised `chapters/15-background.tex`:

- `utherst_decerebrate_2025` — classroom decerebration film (Ben's checkpoint item 5).
  Verified against the video page and full transcript: channel UTHealth Houston School
  of Dentistry, published 2025-04-15, duration 50:13; the cat is exercised first
  (righting, placing reaction), then decerebrated under open-drop ether ("more ether on
  this cat than ... the last horse"), with carotid ligation, tracheal cannula, and
  trephination of the skull. The transection level is also in the transcript
  (auto-caption garbles "colliculi" as "caliculis"): "if we make our transsection
  right here and this is a natural place just between the superior and inferior
  caliculis and actually in the cat following the line of the tentorum cerebelli" —
  i.e., between the superior and inferior colliculi, along the tentorium cerebelli.

```bibtex
@misc{utherst_decerebrate_2025,
	author = {{UTHealth Houston School of Dentistry}},
	title = {1155 {Explicit} {Decerebrate} {Cat} [video]},
	year = {2025},
	howpublished = {\url{https://www.youtube.com/watch?v=1nkczBM5YiM}},
	note = {Classroom demonstration film. Accessed October 2026},
}
```

- `stuart_reflections_2007` — "Reflections on integrative and comparative movement
  neuroscience" (verified: Integrative and Comparative Biology 47(4):482--504,
  doi:10.1093/icb/icm037). Fig. 7 of this paper, per Ben, is the high-decerebrate cat
  treadmill arrangement reprinted from Fig. 1 of Severin et al. 1967; flagged as a
  borrowed-figure candidate for the ethics section.

```bibtex
@article{stuart_reflections_2007,
	author = {Stuart, Douglas G.},
	title = {Reflections on integrative and comparative movement neuroscience},
	journal = {Integrative and Comparative Biology},
	volume = {47},
	number = {4},
	pages = {482--504},
	year = {2007},
	doi = {10.1093/icb/icm037},
}
```

- `severin_1967` — the Russian-language treadmill-decerebrate-cat report whose Fig. 1
  is reprinted as Fig. 7 of Stuart 2007. CAVEAT: exact journal/volume/pages not
  verifiable from here (Ben: "may be in Russian"); confirm against the Stuart 2007
  reference list before finalizing. Stuart 2007 Table 1 cites "Severin et al. 1967"
  alongside Shik et al. 1966 for cat treadmill locomotion. Probable form:

```bibtex
@article{severin_1967,
	author = {Severin, F. V. and Orlovsky, G. N. and Shik, M. L.},
	title = {[Work of the muscles during controlled locomotion]},
	journal = {Biofizika},
	year = {1967},
	note = {In Russian; English translation in Biophysics. Verify exact bibliographic details against Stuart (2007)},
}
```

- `geijtenbeek_scone_2019` — SCONE predictive-simulation software, SimTK-hosted
  (verified: simtk.org/projects/scone; JOSS 4(38):1421).

```bibtex
@article{geijtenbeek_scone_2019,
	author = {Geijtenbeek, Thomas},
	title = {{SCONE}: {Open} {Source} {Software} for {Predictive} {Simulation} of {Biological} {Motion}},
	journal = {Journal of Open Source Software},
	volume = {4},
	number = {38},
	pages = {1421},
	year = {2019},
	doi = {10.21105/joss.01421},
}
```

- `shinohara_mechanisms_2025` — interlimb coordination under sudden ground loss
  (verified via AARL Zotero group, item F55LBFZG). If a journal version appears before
  submission, prefer it over the bioRxiv preprint.

```bibtex
@article{shinohara_mechanisms_2025,
	author = {Shinohara, Kota and Ambe, Yuichi and Kim, Yongi and Mano, Fu and Ramalingasetty, Shravan Tata and Lockhart, Andrew B. and Markin, Sergey N. and Ausborn, Jessica and Rybak, Ilya A. and Danner, Simon M. and Aoi, Shinya},
	title = {Mechanisms of adaptive interlimb coordination to sudden ground loss: a neuromusculoskeletal modeling study},
	journal = {bioRxiv},
	pages = {2025.11.11.687930},
	year = {2025},
	doi = {10.1101/2025.11.11.687930},
}
```

Amendment to an EXISTING entry (Ben's checkpoint item 5; writers cannot edit thesis.bib):

- `decerebrate_cat_footage` (ProofFinal/thesis.bib, line ~1934) — two fixes:
  1. `howpublished = {\\url{...}}` doubles the backslash, printing a literal
     `\url{https://...}` (the "missing space between url and https" symptom). Correct
     to `howpublished = {\url{https://www.youtube.com/watch?v=hAteWL78NMc}}`.
  2. Add the channel as author so the entry cites properly:
     `author = {{WhiteCoatWaste}}` (channel name is one word, WhiteCoatWaste).

## Requested by 50-database.tex (Writer — Database chapter, 2026-10-16)

New keys cited in the revised `chapters/50-database.tex` (justification opening and the
two new tables `tab:sadb_afferents` and `tab:sadb_pathways`). Every entry below was
verified against Crossref by DOI on 2026-10-16, and every paper is a record in the
SADb corpus (SADb_audit knowledge_base pathway/afferent pages carry the DOIs used).
Existing keys reused by the same tables (no action): `dietz_significance_2000`,
`akay_degradation_2014`, `van_de_crommert_review_1998`, `windhorst_muscle_2007`,
`windhorst_spinal_2022`.

- `eccles_integrative_1958` — Ia monosynaptic excitation anchor (Table B; linked to the
  corpus pathway record "Ia monosynaptic excitation", 7 papers).

```bibtex
@article{eccles_integrative_1958,
	author = {Eccles, Rosamond M. and Lundberg, A.},
	title = {Integrative pattern of {Ia} synaptic actions on motoneurones of hip and knee muscles},
	journal = {The Journal of Physiology},
	volume = {144},
	number = {2},
	pages = {271--298},
	year = {1958},
	doi = {10.1113/jphysiol.1958.sp006101},
}
```

- `capaday_amplitude_1986` — phase modulation of the Ia pathway (Table B row 1).

```bibtex
@article{capaday_amplitude_1986,
	author = {Capaday, C. and Stein, R. B.},
	title = {Amplitude modulation of the soleus {H}-reflex in the human during walking and standing},
	journal = {The Journal of Neuroscience},
	volume = {6},
	number = {5},
	pages = {1308--1313},
	year = {1986},
	doi = {10.1523/JNEUROSCI.06-05-01308.1986},
}
```

- `jankowska_effect_1967` — reciprocal organization of Ia pathways (Tables A and B).

```bibtex
@article{jankowska_effect_1967,
	author = {Jankowska, E. and Jukes, M. G. M. and Lund, S. and Lundberg, A.},
	title = {The effect of {DOPA} on the spinal cord {5}. {Reciprocal} organization of pathways transmitting excitatory action to alpha motoneurones of flexors and extensors},
	journal = {Acta Physiologica Scandinavica},
	volume = {70},
	number = {3-4},
	pages = {369--388},
	year = {1967},
	doi = {10.1111/j.1748-1716.1967.tb03636.x},
}
```

- `pratt_ia_1987` — Ia inhibitory interneurons rhythmically active in fictive
  locomotion (Table B row 2).

```bibtex
@article{pratt_ia_1987,
	author = {Pratt, C. A. and Jordan, L. M.},
	title = {{Ia} inhibitory interneurons and {Renshaw} cells as contributors to the spinal mechanisms of fictive locomotion},
	journal = {Journal of Neurophysiology},
	volume = {57},
	number = {1},
	pages = {56--71},
	year = {1987},
	doi = {10.1152/jn.1987.57.1.56},
}
```

- `jankowska_common_1981` — common Ia/Ib interneurones (Table B row 3).

```bibtex
@article{jankowska_common_1981,
	author = {Jankowska, E. and Johannisson, T. and Lipski, J.},
	title = {Common interneurones in reflex pathways from group {1a} and {1b} afferents of ankle extensors in the cat},
	journal = {The Journal of Physiology},
	volume = {310},
	number = {1},
	pages = {381--402},
	year = {1981},
	doi = {10.1113/jphysiol.1981.sp013556},
}
```

- `stephens_short_1996` — non-reciprocal group I inhibition reduced in stance
  (Table B row 3). Crossref volume/pages verified; issue not exposed by the DOI
  record (print is 743(1--2)).

```bibtex
@article{stephens_short_1996,
	author = {Stephens, M. J. and Yang, J. F.},
	title = {Short latency, non-reciprocal group {I} inhibition is reduced during the stance phase of walking in humans},
	journal = {Brain Research},
	volume = {743},
	pages = {24--31},
	year = {1996},
	doi = {10.1016/S0006-8993(96)00977-8},
}
```

- `mccrea_disynaptic_1995` — locomotor reversal of group I action (Table B row 4).

```bibtex
@article{mccrea_disynaptic_1995,
	author = {McCrea, D. A. and Shefchyk, S. J. and Stephens, M. J. and Pearson, K. G.},
	title = {Disynaptic group {I} excitation of synergist ankle extensor motoneurones during fictive locomotion in the cat},
	journal = {The Journal of Physiology},
	volume = {487},
	number = {2},
	pages = {527--539},
	year = {1995},
	doi = {10.1113/jphysiol.1995.sp020897},
}
```

- `pearson_enhancement_1998` — enhancement and resetting of locomotor activity by
  muscle afferents (Table B row 4).

```bibtex
@article{pearson_enhancement_1998,
	author = {Pearson, K. G. and Misiaszek, J. E. and Fouad, K.},
	title = {Enhancement and resetting of locomotor activity by muscle afferents},
	journal = {Annals of the New York Academy of Sciences},
	volume = {860},
	number = {1},
	pages = {203--215},
	year = {1998},
	doi = {10.1111/j.1749-6632.1998.tb09050.x},
}
```

- `dominguez_candidate_2020` — extensor group I reset interneurons inside the
  rhythm-generating layer (Tables A and B, Ib rows).

```bibtex
@article{dominguez_candidate_2020,
	author = {Dom{\'i}nguez-Rodr{\'i}guez, L. E. and Stecina, K. and Garc{\'i}a-Ram{\'i}rez, D. L. and Mena-Avila, E. and Milla-Cruz, J. J. and Mart{\'i}nez-Silva, L. and Zhang, M. and Hultborn, H. and Quevedo, J. N.},
	title = {Candidate interneurons mediating the resetting of the locomotor rhythm by extensor group {I} afferents in the cat},
	journal = {Neuroscience},
	volume = {450},
	pages = {96--112},
	year = {2020},
	doi = {10.1016/j.neuroscience.2020.09.017},
}
```

- `hiebert_contribution_1996` — flexor afferents timing the stance-to-swing
  transition (Table B row 6).

```bibtex
@article{hiebert_contribution_1996,
	author = {Hiebert, G. W. and Whelan, P. J. and Prochazka, A. and Pearson, K. G.},
	title = {Contribution of hind limb flexor muscle afferents to the timing of phase transitions in the cat step cycle},
	journal = {Journal of Neurophysiology},
	volume = {75},
	number = {3},
	pages = {1126--1137},
	year = {1996},
	doi = {10.1152/jn.1996.75.3.1126},
}
```

- `edgley_interneuronal_1987` — midlumbar group I/II relay (Tables A and B).

```bibtex
@article{edgley_interneuronal_1987,
	author = {Edgley, S. A. and Jankowska, E.},
	title = {An interneuronal relay for group {I} and {II} muscle afferents in the midlumbar segments of the cat spinal cord},
	journal = {The Journal of Physiology},
	volume = {389},
	number = {1},
	pages = {647--674},
	year = {1987},
	doi = {10.1113/jphysiol.1987.sp016676},
}
```

- `mazzaro_afferent_2006` — group II stance-phase extensor support in humans
  (Tables A and B).

```bibtex
@article{mazzaro_afferent_2006,
	author = {Mazzaro, N. and Grey, M. J. and Feix do Nascimento, O. and Sinkj{\ae}r, T.},
	title = {Afferent-mediated modulation of the soleus muscle activity during the stance phase of human walking},
	journal = {Experimental Brain Research},
	volume = {173},
	number = {4},
	pages = {713--723},
	year = {2006},
	doi = {10.1007/s00221-006-0451-5},
}
```

- `sherrington_flexion_1910` — the flexion reflex and crossed extension classic
  (Tables A and B, cutaneous rows).

```bibtex
@article{sherrington_flexion_1910,
	author = {Sherrington, C. S.},
	title = {Flexion-reflex of the limb, crossed extension-reflex, and reflex stepping and standing},
	journal = {The Journal of Physiology},
	volume = {40},
	number = {1-2},
	pages = {28--121},
	year = {1910},
	doi = {10.1113/jphysiol.1910.sp001362},
}
```

- `forssberg_stumbling_1979` — the stumbling corrective reaction (Tables A and B).

```bibtex
@article{forssberg_stumbling_1979,
	author = {Forssberg, H.},
	title = {Stumbling corrective reaction: a phase-dependent compensatory reaction during locomotion},
	journal = {Journal of Neurophysiology},
	volume = {42},
	number = {4},
	pages = {936--953},
	year = {1979},
	doi = {10.1152/jn.1979.42.4.936},
}
```

- `prochazka_ensemble_1998` — ensemble afferent coding during locomotion; the signal
  summary behind Table A's Ia/II/Ib rows.

```bibtex
@article{prochazka_ensemble_1998,
	author = {Prochazka, A. and Gorassini, M.},
	title = {Ensemble firing of muscle afferents recorded during normal locomotion in cats},
	journal = {The Journal of Physiology},
	volume = {507},
	number = {1},
	pages = {293--304},
	year = {1998},
	doi = {10.1111/j.1469-7793.1998.293bu.x},
}
```

- `schomburg_flexor_1998` — flexor reflex afferents reset the step cycle (Table A,
  group III/IV row).

```bibtex
@article{schomburg_flexor_1998,
	author = {Schomburg, E. D. and Petersen, N. and Barajon, I. and Hultborn, H.},
	title = {Flexor reflex afferents reset the step cycle during fictive locomotion in the cat},
	journal = {Experimental Brain Research},
	volume = {122},
	number = {3},
	pages = {339--350},
	year = {1998},
	doi = {10.1007/s002210050522},
}
```

## Amendment requested by 15-background.tex (Writer — Background chapter, 2026-10-16 figure round)

- `dirusso_investigating_2023` (ProofFinal/thesis.bib) — DOI correction. The entry's
  `doi = {10.1088/1741-2552/aceb85}` resolves to HTTP 404. The article
  (J. Neural Eng. 20(6):066006, published 16 November 2023, open access CC BY 4.0) is
  published as DOI 10.1088/1741-2552/acfdcc, verified on the IOP article page
  (https://iopscience.iop.org/article/10.1088/1741-2552/acfdcc). Change the entry's
  doi field to:

  doi = {10.1088/1741-2552/acfdcc}

  Reason this surfaced: the chapter now reproduces Figure 5 of this article under
  CC BY 4.0 (sec:bg_sensory), and the caption's license line should point at a
  resolvable DOI.

