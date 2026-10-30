# Publication-note recommendation

I would identify the two SADb presentations clearly, and distinguish the published Morrow–Bolen–Hunt placement work from the new torque manuscript. A presentation note should say **presented as an abstract/poster**, rather than imply that the complete dissertation chapter was published as a journal paper. Merely citing related prior work is different from reproducing or adapting its material.

## PSU requirement and its boundary

PSU requires citation/publication identification when a chapter, major section, or paper incorporates previously published material. For incorporated published co-authored work, its handbook specifies a Chapter Heading Page with full citation information and an author-contribution statement covering all co-authors. It recommends identifying co-authorship/status for unpublished work too. It does not separately declare that every poster appearance makes an entire later chapter a previously published paper. Applying the requirement depends on what material the chapter actually incorporates. [PSU, March 2026 ETD handbook, “Format of the ETD,” “Chapters,” and “Chapter Heading Page & Author Contribution Statement”](https://www.pdx.edu/gradschool/sites/gradschool.web.wdt.pdx.edu/files/2026-03/ETD_Formatting_Requirements_and_Format_Review_Handbook_March_2026.pdf)

## Chapter 5: two verified SADb presentations

The current first paragraph already identifies both meetings at `50-database.tex:4`. I would consolidate that history into one short presentation note, or keep the existing introductory sentence, rather than repeat both. Suggested scope, not an authorized insertion:

> The SADb resource and AI-assisted curation study discussed in this chapter were presented as abstracts and posters at Neuroscience 2023 and Neuroscience 2024, respectively [the two full citations below]. This chapter develops that work and describes the subsequent SADb Explorer application.

Verified records:

| Output | Authors and venue | Existing record and local output |
|---|---|---|
| *Sensory Afferent Database: Advancing towards a comprehensive understanding of sensory afferents and motor control* | B. P. Bolen, K. Deng, and A. J. Hunt; Neuroscience 2023, November 11–15; abstract and poster | `bolen_sensory_2023`; AARL Zotero `Z2ZR9N9B`. Local poster: [2023 SADb poster](</D:/Github/Bipedal_Robot/Documentation/Posters/Bolen et al. - 2023 - Sensory Afferent Database Advancing towards a com.pdf>). |
| *Using AI Tools to Enhance the Sensory Afferent Database: A Comparative Study* | B. P. Bolen and A. J. Hunt; Neuroscience 2024, October 5–9; abstract and poster | `bolen_using_2024`; personal Zotero `LGLRY29Q`, AARL `FSRJIPCN`. Local poster: [SfN24_Bolen.pdf](</D:/Github/Bipedal_Robot/Documentation/Posters/SfN24_Bolen.pdf>). |

An important catalog distinction: personal Zotero `4E92EZ2J` and AARL `7N25LCWZ` describe the **same-title poster at the 2023 NeuroNex Annual Meeting, November 15–16**. They are not the Neuroscience 2023 record. The separate AARL item `Z2ZR9N9B` confirms that SfN presentation, so the current chapter's meeting statement is supported. The current BibTeX entry combines the SfN note with the NeuroNex URL; align the venue/link if revising its metadata. No DOI, abstract number, or additional publication status was invented.

If previously published co-authored text/figures or a substantial section is incorporated, use the PSU citation/contribution treatment for that material. Do not invent contributor roles from author order. Existing presentation authorship is verified; individual roles require their actual record.

## Chapter 4: prior published placement work

The published source is **Connor Morrow, Benjamin Bolen, and Alexander J. Hunt (2020), “Optimization of Artificial Muscle Placements for a Humanoid Bipedal Robot,” in *Biomimetic and Biohybrid Systems*, LNCS 12413, pp. 257–269, DOI 10.1007/978-3-030-64313-3_25**. It is a Living Machines conference-proceedings paper, not the unpublished torque manuscript. The [publisher's proceedings contents](https://link.springer.com/book/10.1007/978-3-030-64313-3?page=2) confirm title, authors, and pages; the existing bibliography key is `morrow_optimization_2020`.

Current Chapter 4 already has an unpublished-manuscript note at `30-torque.tex:5` and cites the earlier method at lines 41, 342, and 344. If the placement subsection materially reproduces/adapts the 2020 method, I would add a narrowly scoped acknowledgment of that source to the existing note and the applicable contribution statement. Do not imply that the new compliance/wrapping identification, physical tests, or later transmission survey appeared in 2020. If the subsection only summarizes and cites the prior method, the citation is appropriate; a separate publication note is not required merely because the paper is related. This is an editorial recommendation, not a claim that the whole chapter is a reprint.

Publication records were found in personal Zotero `5RXI4MTV` and AARL `8CBMCZYH`/`9JFNCZSL`. The AARL records have PDF attachments `UPBEL3MY` and `JLQEE3WZ`; both local files exist and contain 13 pages. [Verified local paper PDF](</C:/Users/Ben Bolen/Zotero/storage/JLQEE3WZ/Morrow et al. - 2020 - Optimization of Artificial Muscle Placements for a.pdf>). Repository source archive: [Morrow optimization paper ZIP](</D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Connor Morrow Optimization Paper/Optimization of Artificial Muscle Attachments for a Biomimetic Humanoid Robot.zip>); poster: [LM202_poster.pptx](</D:/Github/Bipedal_Robot/Documentation/Posters/LM202_poster.pptx>). An archive filename is not substituted for the verified published title.

Both personal and AARL libraries were read through the running Zotero local API. Metadata, actual poster titles/authors, saved-result paths, and source checks are recorded in `TABLE412_PUBLICATION_SOURCE_CHECK.json`. This pass changed no manuscript, bibliography, figure, Zotero item, or scientific result. It created recommendations and a separate ZCode prompt only.
