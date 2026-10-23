from pathlib import Path
import re
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
old=(C/'review/before/50-database.tex').read_text(encoding='utf-8-sig')
def table(kind,label):
    return next(m.group() for m in re.finditer(r'\\begin\{'+kind+r'\}(?:\[[^\]]*\])?.*?\\end\{'+kind+r'\}',old,re.S) if r'\label{'+label+'}' in m.group())
main=r'''\chapter{The Sensory Afferent Database}\label{ch:database}
\graphicspath{{figs/50-database/}}

\section{The Resource}\label{sec:sadb_resource}

The Sensory Afferent Database (SADb) is a curated literature resource for studying sensory feedback in locomotion and posture. It connects papers to afferent classes, experimental preparations, feedback pathways, and computational models. I developed it within the AARL's NeuroNex Communication, Coordination, and Control in Nervous Systems (C3NS) effort and presented the resource at Neuroscience 2023 \citep{bolen_sensory_2023}.

The primary deliverable is a browser application for searching paper records, counting them by selected classifications, and exploring citation neighborhoods (Figure~\ref{fig:sadb_explorer}). Its embedded data support offline use; its separate online-search tab requires a network connection. The distributed snapshot contains 943 audited paper records, including 382 marked as archived after removal from the editable working base. The application distributes paper records and derived indexes, rather than the complete five-table curation base.

The working base in Airtable contains linked Papers, Feedback, Models, Review Papers, and Recorders tables. Papers hold bibliographic information and short mechanism summaries grounded in the source text. Feedback records describe pathways and link them to supporting papers; Models link computational implementations to their source literature. The other tables organize review literature and recording methods. Paper matching uses DOI where available and normalized titles otherwise. These linked records support curation, but the non-paper tables have not yet been exported as complete record sets into the application.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{sadb_explorer.png}
\caption[The SADb Explorer]{SADb Explorer paper-table view. Searches and filters expose bibliographic information, curation notes, preparation tags, and afferent or pathway links. The embedded corpus is available offline.}
\label{fig:sadb_explorer}
\end{figure}

\section{Curation Methods}\label{sec:sadb_methods}

The collection began with my Zotero library, the laboratory's Biology library, and the papers already curated in Airtable. Records were reconciled by DOI and title, then expanded through citation searches and the laboratory's PDF collection. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication.

Inclusion centered on experimental studies, reviews, and computational models of sensory feedback in locomotion and posture. Duplicates and material outside that scope were removed from the working base. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable; missing text was flagged rather than supplied from recollection. The retained and archived sets distinguish the working collection from the broader audit.

\subsection{AI-Assisted Triage}\label{sec:sadb_triage}

I evaluated AI-assisted literature triage in work presented at Neuroscience 2024 \citep{bolen_using_2024}. AI tools helped draft and classify records, but their output required checking against the source. I retained responsibility for inclusion, classification, and deletion decisions. The role of these tools was to assist curation, not to establish the correctness of a physiological claim. Their use does not make the corpus complete or remove the need to inspect the underlying papers.

\section{Coverage and Use}\label{sec:sadb_corpus}

The audited paper set contains 943 records; 561 were retained in the working base. Of those retained papers, 523 have curation notes, 426 have preparation tags, 86 have an afferent classification, and 258 have notes relating their findings to robot or simulation use. These counts describe curation coverage, not the number of independent physiological findings or a measure of their reliability. The more detailed record counts and linked pathway examples are given in Appendix~\ref{app:sadb_records}.

A user can begin with a preparation or afferent class, retrieve the matching papers, and inspect their notes and source links. The pivot view counts papers along two selected dimensions (Figure~\ref{fig:sadb_pivot}); the citation view provides a second route from an anchor paper to related literature (Figure~\ref{fig:sadb_clusters}). These views help locate evidence, while the original experiment determines the conditions under which a finding applies.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.96\textwidth]{sadb_pivot_app.png}
\caption[The SADb Explorer pivot view]{Application pivot view, with afferent classifications as rows and animal preparations as columns. Cells count tagged papers; selecting a cell retrieves its records. Multiple tags allow a paper to contribute to more than one cell.}
\label{fig:sadb_pivot}
\end{figure}

\subsection{Use in Neural-Controller Design}\label{sec:sadb_use}

The resource helped organize the sensory channels used in Chapter~\ref{ch:sns}: spindle-related length and velocity signals, force feedback, and contact events. The linked records provide sources for selecting a pathway and examining its preparation, sign, and phase dependence. They are not ready-made controller parameters. Translating an experimental result into a robot circuit still requires decisions about sensor encoding, species and task differences, and gain selection, followed by verification of the implemented behavior.

The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
'''
(C/'chapters/50-database.tex').write_text(main,encoding='utf-8')
app=r'''\chapter{Appendix E: Supporting SADb Records}\label{app:sadb_records}
\graphicspath{{figs/50-database/}}

This appendix records the curation snapshot and pathway examples supporting Chapter~\ref{ch:database}. Counts refer to the stated audited or retained set; they are not estimates of the size of the neuroscience literature.

\section{Record Sets and Coverage}

Before pruning, the working base contained 1,393 records across its five tables. Removing 436 records left 957: 561 Papers, 39 Feedback, 244 Models, 109 Review Papers, and 4 Recorders. The removed set comprised 382 paper records and 54 linked non-paper records. The application embeds the 943 audited papers, including the archived papers, but does not yet distribute the 450 non-paper records of the full pre-prune base as complete record sets.

Among the 561 retained papers, 523 carry curation notes, 426 preparation tags, 240 notes assessing animal-study potential, 118 links to named feedback pathways, 86 an afferent classification, and 258 robot/simulation translation notes. There are 550 DOIs; 38 records remain without grounding text. These fields overlap and should not be summed as disjoint categories.

The citation graph was constructed over the 943-paper audit: 730 papers were connected by 11,130 directed in-corpus citation pairs, and the Louvain partition produced 18 named groups. Citation proximity supplies a navigation aid, not evidence that two papers establish the same mechanism. Figure~\ref{fig:sadb_clusters} shows the application's paper-level neighborhood view.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.94\textwidth]{sadb_clusters_app.png}
\caption[SADb citation-neighborhood view]{Application citation view focused on Grillner (1985). Soma size represents citation count; axes show year and citation count. Triangles mark 57 in-corpus citers, and circles 13 cited papers. Second-degree papers are gray. The panels identify the selected paper and its connections.}
\label{fig:sadb_clusters}
\end{figure}

\section{Afferent Classes and Linked Pathways}

The tables reproduce the curated signal summaries and representative linked pathways. Their record names are database labels; pathway action depends on the experimental preparation and locomotor state described in the cited sources. Translating these records into a controller requires an explicit encoding and verification step.

'''
t=table('table','tab:sadb_afferents')
t=t.replace(r'\begin{tabularx}{\textwidth}{@{}lllX@{}}',r'\small'+'\n'+r'\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.08\textwidth}>{\raggedright\arraybackslash}p{0.16\textwidth}>{\raggedright\arraybackslash}p{0.22\textwidth}>{\raggedright\arraybackslash}X@{}}')
app+=t+'\n\n'+table('sidewaystable','tab:sadb_pathways')+'\n'
(C/'chapters/96-AppendixE.tex').write_text(app,encoding='utf-8')
p=C/'main.tex';s=p.read_text();s=s.replace(r'\include{chapters/95-AppendixD}',r'\include{chapters/95-AppendixD}'+'\n'+r'\include{chapters/96-AppendixE}');p.write_text(s,encoding='utf-8')
(C/'review/03_SADb.md').write_text('Replaced the chapter with a resource-first description, reproducible curation scope, concise AI role, coverage, and use limits. Full linked-pathway tables, detailed record coverage, and citation-neighborhood screenshot moved to Appendix E. App screenshots replace the retired cluster/pivot PDFs. Distinguishes the app paper export from the five-table working base and the planned full-record export. Removed failed-viewer history, import-batch chronology, general attacks on bibliographies, and the claim that database lookup alone specifies a neural connection. Original counts preserved in main text or Appendix E except omitted operational chronology; audit files retain the complete baseline. Existing image pixels unchanged.\n')
print('Saved SADb chapter and supporting Appendix E')
