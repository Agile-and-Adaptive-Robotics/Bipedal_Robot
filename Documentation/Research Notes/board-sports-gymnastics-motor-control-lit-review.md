# Board Sports & Gymnastics: Biomechanics and Neural Control — Literature Notes

Compiled 2026-09-08 from web-verified academic sources (journal articles, textbooks, peer-reviewed conference proceedings only). Every item in the source lists was confirmed via a real search result or fetched page; secondhand/unverifiable items are explicitly flagged at the end of each section. Companion summary of depth: skateboarding *human biomechanics* is sparse (a 2024 Sports Medicine scoping review found only 8 studies) but its *dynamics/control-modeling* lineage is strong; snowboarding performance biomechanics is sparse-to-moderate with almost no neural-control work; gymnastics is where the human neural-control literature actually lives.

---

## 1. Skateboarding

### Depth assessment
- **Experimental biomechanics: SPARSE.** Diewald et al. (2024, *Sports Medicine*) scoping review located only eight biomechanical studies of skateboarding — almost all lab-based, isolated ollie studies, no ramp/quarter-pipe work.
- **Board–rider dynamics/control modeling: MODERATE (the richest "serious" vein).** Coherent lineage: Hubbard (1979/1980) → Budapest delayed-feedback group (Stépán, Insperger, Várszegi, Molnár) → 2024 optimal-control work (Kogelbauer; Heinen).
- **Neural/motor control specifically: VERY SPARSE.** A few posturography comparisons, one auditory-anticipation EMG study, and the adjacent "human balancing with reaction-time delay" literature. No skateboarding-specific APA or sensory-reweighting studies exist.

### Biomechanics findings
- **Ollie kinetics** (Frederick et al. 2006, *J Appl Biomech*): 7 elite skaters over AMTI force plate; pop forces ~2.25 BW; peak landing vertical GRF 4.52 ± 0.58 BW at 40–50 ms after contact, borne by the forefoot; skaters land deliberately "firm" to stabilize.
- **Landing severity** (via Diewald review, secondhand): handrail landings ~7.98 BW (5344 N); bail-outs up to ~12.09 BW; approach ~4.5 m/s.
- **Ollie kinematics** (Wood et al. 2020, *Apunts*, 11-camera Vicon): no static-vs-rolling difference; greater front-knee flexion ↔ higher board rise. Candotti et al. (2012): ~76% of ollie performance explained by lower-limb power.
- **Modeling/simulation**: Nakashima & Chida (2021, JSME) validated multibody ollie simulation — height needs early rapid front-foot pull-up. Heinen et al. (2024, *Sports Engineering*, TU Delft) direct-collocation optimal control — back-foot-dominant strategy on a smaller board raises ollie height ~12%. Wu et al. (2024, *Front. Bioeng. Biotechnol.*) FE foot/board model — peak forefoot plantar force 830 N, highest stress MT2–MT4.
- **Pumping**: Kogelbauer et al. (2024, *Phys. Rev. Research* 6, 033132) — half-pipe pumping as a variable-length pendulum with friction; optimal crouch/stand timing; analogous to swing pumping.
- **Dynamics/stability**: Hubbard (1979, *J Appl Mech*) — lateral dynamics, higher speed stabilizes roll. Hubbard (1980, *J Biomech*) — rider as rigid body with human control. Kremnev & Kuleshev (2010, DCDS-S) — nonlinear analysis, no rider. Rosatello et al. (2015, ASME IDETC) — wobble onset depends predominantly on HUMAN CONTROL characteristics, not board geometry. **Várszegi et al. (2016, *J R Soc Interface*) — speed wobble emerges from reflex delay in the human control loop** (key result for anyone modeling delayed-feedback balance).

### Motor control findings
- Turkish dynamic-posturography study (institutional record confirmed, authors not retrievable): skateboarders better overall Sensory Organization Test balance (p<0.001), differences in visual and vestibular conditions.
- Cesari et al. (2014, *PLoS ONE*): board/landing SOUND affects action anticipation and muscle activation in a simulated ollie — acoustic cues time motor output.
- Expert–novice: Walsh, Creekmur & Wojcik (ISBS 2006) — beginning vs skilled ollie force–time profiles. Vorlíček et al. (2015, *Gymnica*) — switch-stance ollies raise preparatory muscle activity in back limb. Thompson et al. (2025, *Applied Sciences*) — vastus lateralis dominant activation. Ab Rasid et al. (2024, *PLoS ONE*) — stork/star-excursion balance + strength best discriminate skill (90% classification). Vargas et al. (2015) — stance direction (regular/goofy) shapes directional balance.
- **Adjacent, rich**: Molnár, Zelei & Insperger (2021 *J R Soc Interface*; 2022 *J Biomech*) — rolling balance board as skateboard-roll-plane task to estimate human reaction-time delay; "critical delay" as task-difficulty measure for frontal-plane balancing.

### Robotics crossover
- **Chen, Rogers, Zhang & Sreenath (2019), arXiv:1907.11353 — "Feedback Control for Autonomous Riding of Hovershoes by a Cassie Bipedal Robot"**: trajectory-optimization-based feedback; Cassie balances on two hoverboards, regulates forward/yaw velocity, turns, traverses slopes/stairs/rough terrain; visual SLAM + Dijkstra planning. No other confirmed skateboard-riding robot papers.

### Sources (confirmed)
1. Hubbard (1979), *J Appl Mech* 46(4):931–936. https://asmedigitalcollection.asme.org/appliedmechanics/article/46/4/931/422552/
2. Hubbard (1980), *J Biomech* 13:745–754. https://pubmed.ncbi.nlm.nih.gov/7440589/
3. Kremnev & Kuleshev (2010), *DCDS-S* 3(1):85. https://www.aimsciences.org/article/doi/10.3934/dcdss.2010.3.85
4. Rosatello et al. (2015), ASME IDETC/CIE. https://hal.science/hal-01369978
5. Várszegi, Takács, Stépán, Hogan (2016), *J R Soc Interface* 13(121):20160345. https://doi.org/10.1098/rsif.2016.0345
6. Kogelbauer, Koyama, Callan, Shinomoto (2024), *Phys. Rev. Research* 6:033132. https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.033132
7. Frederick, Determan, Whittlesey, Hamill (2006), *J Appl Biomech* 22(1):33–40. https://doi.org/10.1123/jab.22.1.33
8. Wood, Oliveira, Santos, Rodacki, Lara (2020), *Apunts* 141:87–91. https://revista-apunts.com/en/3d-kinematic-analysis-of-the-ollie-maneuver-on-the-skateboard/
9. Nakashima & Chida (2021), *Mech Eng J (JSME)* 8(5):21-00230. https://www.jstage.jst.go.jp/article/mej/8/5/8_21-00230/_article
10. Heinen et al. (2024), *Sports Engineering* 27. https://doi.org/10.1007/s12283-023-00448-y
11. Wu, Wang, Deng, Guo, Zhu (2024), *Front Bioeng Biotechnol* 12:1382161. https://doi.org/10.3389/fbioe.2024.1382161
12. Diewald, Neville, Cronin, Read, Cross (2024), *Sports Medicine* 54(6):1399–1418 (scoping review). https://pmc.ncbi.nlm.nih.gov/articles/PMC11239769/
13. Walsh, Creekmur, Wojcik (2006), ISBS Proc. 24. https://ojs.ub.uni-konstanz.de/cpa/article/view/265
14. Cesari, Camponogara, Papetti, Rocchesso, Fontana (2014), *PLoS ONE* 9(3):e90156. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090156
15. Vorlíček et al. (2015), *Gymnica*. https://gymnica.upol.cz/pdfs/gym/2015/01/06.pdf
16. Thompson et al. (2025), *Applied Sciences* 15(15):8280. https://www.mdpi.com/2076-3417/15/15/8280
17. Ab Rasid et al. (2024), *PLoS ONE* 19(2):e0296467. https://pmc.ncbi.nlm.nih.gov/articles/PMC10852284/
18. Molnár, Zelei, Insperger (2021), *J R Soc Interface* 18(176):20200956. https://doi.org/10.1098/rsif.2020.0956
19. Molnár, Zelei, Insperger (2022), *J Biomech*. https://www.sciencedirect.com/science/article/pii/S0021929022001683
20. Vargas, Wong, Patton, Brown (2015), conference (CSU Fullerton) — via ResearchGate only, weak verification. https://www.researchgate.net/publication/297918854
21. Chen, Rogers, Zhang, Sreenath (2019), arXiv:1907.11353. https://arxiv.org/abs/1907.11353
22. Candotti et al. (2012), via SciELO. https://www.scielo.br/j/rbce/a/dbKJJmKFdyYjcSNfjcpnyXb/

**Flagged secondhand (cited inside Diewald review, not independently verified):** Determan et al. 2010; Leuchanka et al. 2017; Ou et al. 2021. Authors unretrievable (blocked pages): Medipol posturography study; Kinesiologia Slovenica sway-test paper; full author list of Thompson et al. 2025.

---

## 2. Snowboarding

### Depth assessment
- **Performance biomechanics: SPARSE-TO-MODERATE.** ~15–25 genuine studies in three clusters: (a) 2012–2014 Austrian/Dutch inverse-dynamics wave (Klous; Krüger), (b) active 2024–2026 German/Swiss wearables program (Merz, Gorges, Thelen, Kersting) in Sensors/JSS/ISBS, (c) scattered Polish/Korean balance work. Sample sizes n = 1–12. Nothing comparable to the alpine-skiing literature.
- **Neural/motor control: SPARSE.** One true neurophysiological study (EEG+EMG, Ouyang 2023); rest is force-plate posturography with mostly null results. Expert–novice comparisons nearly absent.

### Biomechanics findings
- **Carved-turn joint loading** (Klous, Müller & Schwameder 2014, *Comput Math Methods Med*): mobile Kistler plates under binding parts + panning-tilting-zooming cameras + inverse dynamics (extended Yeadon model). Steering leg: ankle plantarflexion moments 2.5–3.6 N·m/kg, knee flexion-extension moments 2.3–3.3 N·m/kg (peaks to ~6), ankle longitudinal forces 0.4→0.8 BW through the turn. **Joint moments consistently greater in snowboarding than skiing.**
- **Measurement rig** (Krüger et al. 2012, *Proc IMechE H*): custom two-force-plate + full-body IMU + OpenSim multi-segment inverse dynamics; validated only for intra-individual technique/equipment comparison.
- **Carve physics**: Wu, Igci, Andreopoulos & Weinbaum (2006, *MSSE*) — lift mechanics on compressible snow; Wu & Sun (2011, *MSSE*) — board-width function, multiple rider equilibrium positions, pore-pressure distributions (the main board-optimization theory strand).
- **Halfpipe/big air** (German cluster): Thelen et al. (ISBS 2024) GNSS+IMU+pressure insoles, 1 rider/32 airs — COM velocity 5.0–13.0 m/s, back-leg knee flexion 29–64°, peak rear-foot load 1269 ± 123 N. Merz et al. (2025 JSS; 2026 *Int J Sports Sci Coaching*) — deterministic models of airtime vs angular-velocity trade-off ("go big or spin fast"), Beijing-2022 halfpipe. Jiang et al. (2025, *Appl Sci*) — takeoff kinematics/kinetics differentiating trick difficulty. Gorges et al. (2024, *Sensors*) — U-Net on IMU beats threshold algorithms for airtime detection.
- **Terrain park**: Kurpiers et al. (2017, *Injury*) — video analysis of 704 jumps: fall once per 5 jumps, predicted by flat/knuckle landings, poor board orientation, spinning tricks.
- **Other measurement**: ISBS COM-trajectory estimation in giant slalom (2020, art. 158); TOF-camera 3D vision capture 98.6% (Li et al. 2021); IMU ankle rig on lab simulator (Park et al. 2015); grab detection via jerk (ISBS 2024 art. 173).

### Motor control findings
- **The one neurophysiology paper** — Ouyang & Chen (2023, *Percept Mot Skills*): 12 trained halfpipe riders vs 12 untrained, 30/60 cm drops. Pre-drop EEG: increased frontal theta, reduced central/parietal alpha with height, larger in athletes — interpreted as stronger **predictive (feedforward) sensorimotor modulation**. Landings: trained riders higher rectus femoris/gastrocnemius activation with REDUCED antagonist coactivation (higher GAS/TA ratio at 60 cm) — sport-specific neuromuscular patterning.
- Posturography mostly null: Kłos et al. (2019, *Acta Bioeng Biomech*) — no COP-sway change after 6-day course; Staniszewski et al. (2017, *Pol J Sport Tourism*) — no stability gain after 9 days. Staniszewski et al. (2016, *Human Movement*) — stance laterality ↔ lateral preference. Jeon & Eom (2021) — physique/fitness predict national-team balance. Vernillo et al. (2016/2017) — front/rear-leg strength asymmetry, vastus lateralis architecture changes in elites.
- **No snowboard-specific sensory-contribution or feedforward-vs-feedback studies — that subfield is empty.**

### Robotics crossover
- 2003 JSME Hokuriku-Shinetsu conference paper "3-dimensional simulation of snowboarding and development of a snowboard robot" (authors unverifiable, ResearchGate only). Hament & Cater (2017, IEEE) — VR + motion-platform training simulator. Adjacent ski-robot control literature is transferable: multilayer stability control on unknown slopes (IEEE 6094546); carving-turn control (IEEE 5354807); DRL biped ski robot (IEEE 9539926); humanoid ski robot carving via torso-angle control with gate vision (*Sensors* 22(3):816). Zero arXiv "snowboard + robot" papers.

### Sources (confirmed)
1. Klous, Müller, Schwameder (2014), *Comput Math Methods Med* 2014:340272. https://pmc.ncbi.nlm.nih.gov/articles/PMC4181787/
2. Krüger, McAlpine, Borrani, Edelmann-Nusser (2012), *Proc IMechE H* 226(2):170–175. DOI 10.1177/0954411911426938
3. Wu, Igci, Andreopoulos, Weinbaum (2006), *Med Sci Sports Exerc*. https://pubmed.ncbi.nlm.nih.gov/?term=Wu+lift+mechanics+downhill+skiing
4. Wu, Sun (2011), *Med Sci Sports Exerc* 43(10):1955–1963 (via Europe PMC).
5. Thelen, Merz, Gorges, Goldmann, Donath, Kersting (2024), ISBS 42(1):87. https://commons.nmu.edu/isbs/vol42/iss1/87/
6. Merz, Naundorf, Gorges, Kersting (2026), *Int J Sports Sci Coaching*. https://journals.sagepub.com/doi/abs/10.1177/17479541261439267
7. Merz et al. (2025), *J Sports Sci* 43(3):299–307. https://www.tandfonline.com/doi/abs/10.1080/02640414.2025.2453788
8. Jiang, Chen, Gao, et al. (2025), *Appl Sci* 15(12):6618. https://www.mdpi.com/2076-3417/15/12/6618
9. Gorges, Davidson, Boeschen, Hotho, Merz (2024), *Sensors* 24(21):6773. https://www.mdpi.com/1424-8220/24/21/6773
10. Gorges & Thelen & Merz (2024), ISBS 42(1):160; Friedl, Gorges, Merz (2024), ISBS 42(1):173. https://commons.nmu.edu/isbs/vol42/iss1/160/
11. ISBS (2020) 38(1):158 — COM trajectory, snowboard GS. https://commons.nmu.edu/isbs/vol38/iss1/158/
12. ISBS (2017) 35(1):42 — snowboard cross start. https://commons.nmu.edu/isbs/vol35/iss1/42/
13. Ouyang, Chen (2023), *Percept Mot Skills* 130(2):844–862. DOI 10.1177/00315125221148637
14. Kłos, Giemza, Dziuba-Słonina (2019), *Acta Bioeng Biomech* 21(1):97–101.
15. Staniszewski, Zybko, Wiszomirska (2016), *Human Movement* 17(2):119–125. DOI 10.1515/humo-2016-0015
16. Staniszewski, Zybko, Wiszomirska (2017), *Pol J Sport Tourism* 24(2):97–101. DOI 10.1515/pjst-2017-0010
17. Jeon, Eom (2021), *J Exerc Sci Fit* (via PubMed).
18. Falda-Buscaiot, Hintzy (2015), *Comput Methods Biomech Biomed Engin* 18(sup1):1936–1937 (abstract). DOI 10.1080/10255842.2015.1069576
19. Kurpiers, McAlpine, Kersting (2017), *Injury* 48(11):2457–2460. DOI 10.1016/j.injury.2017.08.052
20. Vernillo, Pisoni, Thiébat (2016), *Clin J Sport Med* (via PubMed).
21. Olivié et al. (2025), *Movement & Sport Sciences* — knee kinematics backside vs frontside. https://mbj.episciences.org/en/articles/14562/download
22. Li, Wang, Zhang, Zhang (2021), Wiley art. 8517771 — TOF 3D vision capture. https://onlinelibrary.wiley.com/doi/10.1155/2021/8517771
23. Hament, Cater, et al. (2017), IEEE. https://ieeexplore.ieee.org/document/7992668/
24. Park et al. (2015), *J Biomed Eng Research* (Korea). https://koreascience.kr/article/JAKO201530861235324.image

**Not found / excluded:** no Sports Engineering snowboard-turn EMG paper exists; 2003 JSME snowboard-robot authors unverifiable; "Back et al. 2013" EMG carving paper exists but in a low-tier venue (excluded).

---

## 3. Gymnastics (fallback — where the neural-control literature actually is)

### Depth assessment
- **Neural/motor control: MODERATE, approaching rich.** Dedicated author clusters (Asseman/Caron/Crémieux; Vuillerme & Nougier; Busquets/Federolf; Heinen; Natrup/Wagner; von Laßberg), a Sports Medicine review entry, and a 2025 scoping review. Caveats: no gymnast-specific APA literature; vestibular-adaptation evidence thin and partly negative.
- **Biomechanics generally: RICH.** Dedicated journal (*Science of Gymnastics Journal*, Univ. Ljubljana), IOC handbook chapter, decades of ISBS proceedings, mature somersault/aerial-mechanics line (Yeadon).

### Neural/motor control findings (priority)
- **Superior postural control is real but task-specific**: Asseman et al. (2004 *Neurosci Lett*; 2008 *Gait Posture*) — elite gymnasts' advantage shows mainly in gymnastics-specific stances (beam-like/narrowed) and under reduced/modified visual & proprioceptive inputs; does not fully transfer to ordinary postures. Asseman et al. (2005, *IJSM*) — vision removal degrades gymnasts' sway as task difficulty rises.
- **Automaticity**: Vuillerme & Nougier (2004, *Brain Res Bull*) — dual-task paradigm: gymnasts need fewer attentional resources for sway regulation. Isableu et al. (2017, *Front Hum Neurosci*) — no difference in conventional sway metrics vs other athletes, but significantly higher COP sample entropy in gymnasts = more automatic, less attention-invested control that conventional posturography misses.
- **Sensory reweighting (direct evidence)**: Busquets et al. (2018, *Gait Posture*) — sensory-reweighting capability for posture modulated by both age and gymnastic expertise. Busquets et al. (2021, *Front Psychol*) — gymnastics experience enhances multi-segmental coordination during PROPRIOCEPTIVE reweighting (PCA/Federolf method).
- **Inverted stance**: Gautier, Thouvarecq & Chollet (2007, *J Sports Sci*) — visual & postural control of handstand in experts; 2009 follow-up (*Hum Mov Sci*) — coordination/neuromuscular strategy changes with expertise. Asseman & Gahéry (2005, *Neurosci Lett*) — head position and vision effects on handstand balance.
- **Vestibular adaptation — contested**: von Laßberg, Campos & Beykirch (2020, *PLoS ONE*, 3-year longitudinal) — NO systematic VOR adaptation in elite gymnasts; eye movements in real twisting somersaults (SCEMs) likely proprioceptively/centrally driven. van der Veen et al. (2022, *Front Sports Act Living*) — skill-related adaptive modifications of gaze stabilization (greater VOR suppression) with experience. Net: active, task-embedded vestibulo-ocular strategies adapt; passive reflex gain largely does not.
- **Spotting / aerial gaze**: Heinen (2011, *Int J Sport Psych*) — experimental evidence gymnasts visually spot in back aerial somersaults even with short flight times. Natrup et al. (2020, *Hum Mov Sci*) — trampoline back tucks, higher performers fixate the bed differently. Natrup et al. (2021, *Hum Mov Sci*) — spotting with fixations even during full-twisting somersaults.
- **Feedforward pre-landing activation & stiffness tuning**: Janshen (2000, ISBS) — all muscles pre-activate before touchdown; longer pre-activation in highly trained gymnasts → increased knee/ankle stiffness and redistributed foot pressures. Niespodziński et al. (2021, *J Hum Kinet*) — gymnasts show LOWER normalized pre-landing EMG (≈half the multifidus/gastrocnemius of untrained) yet 13% higher peak GRF and faster force development — tuned, economical co-activation. Pavlasová et al. (2025, *Front Sports Act Living*, scoping review, 8 studies) — pre-landing onset 80–90 ms (vastus lateralis) before contact; adults (not children) modulate whole-body stiffness to the surface; "soft" (>63° knee flexion) vs "stiff" (<63°) taxonomy; GRFs 7.1–15.8 BW.
- **"Zeroing"**: Pain et al. (2007) — landing goal framed as reducing mass-center velocity to zero without a step/hop/jump (Portsmouth portal abstract).

### Biomechanics findings
- Beam: back handspring angular-momentum generation/control on beam (*J Biomech*, 25 gymnasts; title/venue confirmed, authors not retrieved).
- Vault/springboard: Greenwood (1996, ISBS) — direct force measurement of vault takeoff (Kistler); Atiković (2012, *Sci Gymnastics J*) — vault-flight regression models; a 2026 *Scientific Reports* tucked-front-somersault study includes bilateral activation asymmetries.
- **Aerial mechanics (Yeadon cluster)**: Yeadon & Hiley (2014, *J Biomech*) — twisting-somersault control via simulation: twist initiated by contact torque at takeoff OR aerially via asymmetrical arm/hip motion; mid-flight corrections via symmetrical arm/hip adjustments. Companion "Biomechanics of Twisting Somersaults Part II: Contact Twist" (Loughborough repository; venue unconfirmed).
- Landing surfaces: McNitt-Gray, Yokoi & Millward (1994, *J Appl Biomech*) — landing strategies on different surfaces.

### Sources (confirmed)
1. Asseman, Caron, Crémieux (2004), *Neurosci Lett* 358(2):83–86. https://pubmed.ncbi.nlm.nih.gov/15026154/
2. Asseman, Caron, Crémieux (2008), *Gait Posture* 27(1):76–81. https://pubmed.ncbi.nlm.nih.gov/17337190/
3. Asseman, Caron, Crémieux (2005), *Int J Sports Med* 26(2):116–119. https://pubmed.ncbi.nlm.nih.gov/15726495/
4. Asseman, Gahéry (2005), *Neurosci Lett* 375(2):134–137. https://pubmed.ncbi.nlm.nih.gov/15670656/
5. Vuillerme, Nougier (2004), *Brain Res Bull* 63:161–165. https://pubmed.ncbi.nlm.nih.gov/15130706/
6. Isableu, Hlavackova, Diot, Vuillerme (2017), *Front Hum Neurosci* 11:317. https://www.frontiersin.org/journals/human-neuroscience/articles/10.3389/fnhum.2017.00317/full
7. Hrysomallis (2011), *Sports Medicine*. https://pubmed.ncbi.nlm.nih.gov/21395364/
8. Busquets et al. (2018), *Gait Posture* 63:177–183. https://pubmed.ncbi.nlm.nih.gov/29763813/
9. Busquets, Ferrer-Uris, Angulo-Barroso, Federolf (2021), *Front Psychol* 12:661312. https://pubmed.ncbi.nlm.nih.gov/33935920/
10. von Laßberg, Campos, Beykirch (2020), *PLoS ONE*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7735588/
11. van der Veen et al. (2022), *Front Sports Act Living*. https://www.frontiersin.org/journals/sports-and-active-living/articles/10.3389/fspor.2022.824990/full
12. Heinen (2011), *Int J Sport Psychology*. https://pubmed.ncbi.nlm.nih.gov/21628729/
13. Natrup et al. (2020), *Hum Mov Sci*. https://pubmed.ncbi.nlm.nih.gov/32217208/
14. Natrup et al. (2021), *Hum Mov Sci*. https://www.sciencedirect.com/science/article/abs/pii/S0167945720306047
15. Gautier, Thouvarecq, Chollet (2007), *J Sports Sci*. https://pubmed.ncbi.nlm.nih.gov/17654239/
16. Janshen (2000), ISBS Proc. https://ojs.ub.uni-konstanz.de/cpa/article/view/2300
17. Niespodziński et al. (2021), *J Hum Kinet* 78:15–28. https://pmc.ncbi.nlm.nih.gov/articles/PMC8120959/
18. Pavlasová, Bizovská, Gába, Farana, Janura (2025), *Front Sports Act Living*. https://pmc.ncbi.nlm.nih.gov/articles/PMC12179790/
19. McNitt-Gray, Yokoi, Millward (1994), *J Appl Biomech* 10(3):237–252. https://journals.humankinetics.com/view/journals/jab/10/3/article-p237.xml
20. Yeadon, Hiley (2014), *J Biomech* 47(6):1340–1347. https://www.sciencedirect.com/science/article/abs/pii/S0021929014000979
21. Yeadon, "Biomechanics of Twisting Somersaults Part II: Contact Twist" (title confirmed, venue unconfirmed). https://repository.lboro.ac.uk/articles/journal_contribution/The_biomechanics_of_twisting_somersaults_Part_II_contact_twist/9624797
22. **Textbooks**: Winter DA (2009) *Biomechanics and Motor Control of Human Movement*, 4th ed., Wiley; Shumway-Cook A & Woollacott MH (2007) *Motor Control: Theory and Practical Applications*, 3rd ed., LWW; Zatsiorsky VM, ed. (2000) *Biomechanics in Sport* (IOC Medical Commission), incl. Brüggemann's "Biomechanics of gymnastics" chapter — free PDF: https://stillmed.olympics.com/media/Document%20Library/OlympicOrg/IOC/Who-We-Are/Commissions/Medical-and-Scientific-Commission/Encyclopaedia/2000_Zatsiorsky.pdf

**Flagged/unconfirmed:** DeVita & Skelly 1992 venue; beam back-handspring paper authors; no gymnast-specific APA studies found.

---

## 4. Cross-cutting themes & relevance to bipedal-robot / SNS work

1. **Delayed feedback as the destabilizer.** The strongest quantitative control result in board sports is Várszegi et al. 2016: skateboard speed wobble emerges from reflex delay in the human control loop (not board geometry — Rosatello 2015). The Budapest group's balance-board "critical delay" (Molnár & Insperger 2021/2022) turns this into a task-difficulty metric. Directly maps onto reflex-loop design in a synthetic nervous system: loop latency bounds the achievable balance bandwidth.
2. **Feedforward stiffness tuning beats brute-force co-contraction.** Expert landers (gymnastics: Niespodziński 2021, Janshen 2000; snowboarding: Ouyang 2023) use LOWER overall pre-activation with better-targeted activation and reduced antagonist coactivation while absorbing higher loads — economical, predictive stiffening, not maximal co-contraction.
3. **Sensory reweighting, not sensory superiority.** Gymnast expertise shows up as flexible reweighting between visual/proprioceptive/vestibular channels (Busquets 2018/2021) and more automatic control (Isableu 2017), while passive vestibular reflex gain does NOT adapt (von Laßberg 2020). Active gaze strategies (spotting; VOR suppression) do adapt.
4. **Task specificity of expertise.** Balance skill is posture-specific (Asseman 2004/2008) — relevant when choosing a validation task for a controller.
5. **Robotics crossover is thin but pointed.** Cassie on hovershoes (Chen/Sreenath 2019) is the one serious biped-on-board result: TO-based feedback control. Human results suggest a control architecture split: slow feedforward trajectory (line choice, pump timing — Kogelbauer 2024 pendulum-pumping) + fast delayed-feedback roll stabilization with latency-aware limits.
