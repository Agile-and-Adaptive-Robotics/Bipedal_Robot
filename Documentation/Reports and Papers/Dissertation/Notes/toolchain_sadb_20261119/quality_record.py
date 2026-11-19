"""Save Ben's rulings, check exact revised passages, and extend the single reader review table."""
from pathlib import Path
from zipfile import ZipFile
import json, html, re, sys
from bs4 import BeautifulSoup
here=Path(__file__).resolve().parent
diss=here.parents[1]
chap=diss/"Final/chapters"
root=diss.parents[2]
history=here.parent/"advisor_feedback_history"
sys.path.insert(0,str(history))
import review

rulings=[
("BOLEN-NOV19-BIW2LVAR","I told you many times to stop calling one of the neural models \"Hunt\".","Name Ben's family BiW2Lvar+. Give configuration, cell count, and experiment when results are from different versions. Attribute the four-stage calibration method to Hunt and colleagues.","contribution_and_credit"),
("BOLEN-NOV19-DRAWING","Chapter 6 should be the SNS-Simscape block programming.","Show which Simulink blocks are wired and built in Chapter 6. Give the optional connectome drawing one brief appendix note. Put a proposed drawing-to-simulator application in Chapter 7 future work.","reader_use"),
("BOLEN-NOV19-SUSPENDED","That makes it sound like bicycling.","Describe stepping-like movements in the suspended model. Describe what happened when the model was lowered and support was removed as a separate observation.","definitions_and_terminology"),
("BOLEN-NOV19-FACTORS","I wanted you to remove \"Isometric\" because we talk about dynamic factors in future work.","Use Joint Design Factors for the section spanning fixed-knee experiments and proposed dynamic power measurements.","scope_and_uncertainty"),
("BOLEN-NOV19-WORKFLOW","the substance would be changed because figures made for slides aren't always what figures in manuscripts should be","Rebuild the workflow for a reader choosing an experiment. Include the intervention, fixed conditions, signals, and decision. Measure arrows and labels at print size.","evidence_and_validation")]
path=history/"EDITORIAL_EXAMPLES.json"
content=path.read_text(encoding="utf-8")
rows=json.loads(content)["examples"]
new=[]
for key,quote,lesson,tag in rulings:
    if key not in {r["id"] for r in rows}:
        new.append(dict(id=key,source=(here/"REVIEW_TABLE.html").relative_to(root).as_posix(),anchor="Ben's November 19 manuscript annotations and corrections",source_text=quote,application=lesson,status="Ben's correction. No advisor event or deduction.",tags=[tag]))
if new:
    prefix='  "examples": [\n'
    rendered=",\n".join("\n".join("    "+s for s in json.dumps(r,ensure_ascii=False,indent=2).splitlines()) for r in new)
    content=content.replace(prefix,prefix+rendered+",\n",1)
    assert len(json.loads(content)["examples"])==len(rows)+len(new)
    path.write_text(content,encoding="utf8")

# Each locator applies to stable working TeX, including the earlier wording captured in the live archive.
findings=[
("Ben and database reviewer","50-database","Choosing a sensory connection",4,"Why create SADb?","Names the nerve, response, preparation, and modeling choice that make a source-linked collection necessary.","gap_and_motivation"),
("Reader reproducing the lookup","50-database","Open SADb Explorer in a browser",10,"Which app can I open?","The tracked expanded app contains all 1,393 records and exposes the record and archive filters. The dist version only exposes 943 Papers.","reader_use"),
("Database reviewer","50-database","Counting both citing and cited",49,"Which papers belong to the citation graph?","730 is the count of papers with outgoing citations. The union of citing and cited endpoints is 816. The edge count of 11,130 and the 18 communities are unchanged.","evidence_and_validation"),
("Ben and independent reader","50-database","The Table tab retrieves both",51,"Can I follow a pathway record back to a study?","The actual Table search returns six Papers. The item for McCrea et al. links to its DOI and displays Cat, Ia, Ib, and both pathway names. Enlarged genuine screen regions replace the Deng screenshot.","reader_use"),
("Primary experiment and independent reader","50-database","In decerebrate cats at rest",69,"What can a modeler take from the experiment?","The same group I nerve stimulus inhibited the synergist motor pool at rest and excited it during extension in midbrain-evoked fictive locomotion. A signed diagram conveys the net response without assigning unmeasured interneurons or gains.","mechanism"),
("Ben and database reviewer","50-database","Pivot answers a different question",72,"Will a reader understand the app?","The Ib and Cat cell counts 21 Papers. Selecting it opens the contributing rows. The screenshot also shows the return control. The worked Table, Pivot, and Bubble map tasks follow one another.","reader_use"),
("Database reviewer","50-database","two-link display contains 555",64,"What does the map count and which control changes it?","B changes depth and E marks the resulting neighborhood. The 555 displayed Papers include the focus record and 554 other Papers.","evidence_and_validation"),
("Ben and technical reader","60-sns","In the SNS-Simscape interface",112,"Which graphical tool changes a runnable circuit?","SNSGui builds a nonspiking Simulink circuit. The knee and arm builders add the mechanics. The separate drawing tool is absent from Chapter 6.","contribution_and_credit"),
("Defense-folder audit","60-sns","a separate bounded least-squares fit",175,"What did the fitting tool vary?","TPE proposes rhythm and timing parameters. Bounded least squares fits pattern-formation to motor-neuron weights in s3k and syn6.","mechanism"),
("Defense-folder audit","60-sns","a recorder that saves",200,"What can I recover from a completed trial?","Saved timestamps align neural and mechanical signals. Reloading the saved poses reproduces the observed motion.","reader_use"),
("Ben and source profiles","60-sns","The second timing test used",259,"Whose model, and which configuration?","BiW2Lvar+ identifies Ben's controller family. The timed configurations are identified by neuron count, balance feedback, and Renshaw inhibition. The configuration for the reported suspended test is distinct.","contribution_and_credit"),
("Ben and Overleaf author capture","70-overall","The legs performed stepping-like",61,"What happened in suspension and on release?","Retains Ben's observed stepping-like motion in air and fall on removal of support. Removed the repeated conjunction from the live edit.","definitions_and_terminology"),
("Ben and archive","70-overall","Exploded longitudinal half-section",35,"Where is the original FEA illustration?","Restored the original exploded view as Figure 7.2C. Panels A and B retain the axial profile and experimental load cases. The section describes a proposed finite-element investigation.","reader_context"),
("Ben and curation evaluation","70-overall","Future researchers can revisit large language models",57,"What has to change for assisted classification to be viable?","Propose source-grounded labels, an independent curation comparison, and a human decision. The Chomsky, Roberts, and Watumull guest essay is an argument about explanation. The empirical test uses the project's screening experience.","scope_and_uncertainty"),
("Ben and technical reader","70-overall","Future researchers can develop an application",67,"Can the current drawing change what runs?","A diagram-to-simulator application is future work. A direct change in the controller gives a reference for testing the proposed application.","reader_use"),
("Ben and defense source","70-overall","locate the first change in coordination",0,"Which comparison does the slide suggest for a manuscript?","The new portrait Figure 7.5 takes one baseline into two alternative controlled repeats and identifies the signals to compare. Chapter 7.3 supplies the interpretation.","evidence_and_validation"),
("Ben and technical reader","95-AppendixD","A saved connectome drawing",457,"Which controller is drawn?","One note identifies an earlier BiW2Lvar+ configuration and its four-cell rectus femoris and hamstring circuit on each side. Editing the drawing changes the diagram.","definitions_and_terminology"),
("Ben and source model","95-AppendixD","For the recorded Gait2392 body trial",457,"Do the plots and replay correspond to one trial?","The compiled model and saved positions reproduce a view of that trial, with neural, joint, foot-load, and harness signals on their recorded time bases.","reader_use"),
("Ben and independent reader","96-AppendixE","Of the 561 retained Papers",12,"What does Appendix E add?","Distinguishes 957 maintained from 436 archived, and overlapping classifications from unique records. The 426 records have animal or human tags, not preparation tags. Tables E1 and E2 give receptors, possible inputs, and eight named pathways.","evidence_and_validation"),
("Primary experiment and independent reader","96-AppendixE","In the acute spinal preparation",20,"What alternative condition matters to a specialist?","In the clonidine-treated acute spinal preparation, inhibition diminished and the short-latency excitation was absent. A model of that preparation requires a different test.","scope_and_uncertainty")]

observed=review.signals()
index=review.EvidenceIndex(observed)
chapters=[chap/(c+".tex") for c in sorted({f[1] for f in findings})]
base=ZipFile(here/"LIVE_CHAPTER_BASELINE.zip")
cases=[]
for i,(who,c,needle,n,question,answer,tag) in enumerate(findings,1):
    units=[v for v in review.passages(chap/(c+".tex")) if needle in v[-1]]
    assert len(units)==1,(needle,len(units))
    line,_,kind,source=units[0]
    earlier=base.read(c+".tex").decode().splitlines()
    old=earlier[n-1] if n else "No equivalent Chapter 7 figure in the captured manuscript."
    cases.append(dict(id=f"TOOL-SADB-{i:02d}",file=c+".tex",line=line,author=who,question=question,before=old,after=source,after_sha256=review.digest(source),judgment=answer,tags=[tag],historical_examples=index.nearest(source),local_readiness=3))
record=dict(advisor_event_ids_at_review=[r["id"] for r in observed],rubric="GPT-6 Codex editorial review. Read historical comments and current contexts, identify the mechanism, actor, intervention, and warranted inference, then check against primary source and rendered page. Local readiness 3 means source and context checked. 4 requires author or advisor confirmation. Unmarked passages are unlabeled.",passages=cases)
for name,data in [("PASSAGE_REVIEW.json",record)]:
    (here/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
out=review.make_review(chapters,observed,record)
assert not out["saved_judgments_stale_or_absent"] and not out["reopened_judgments"]
(root/"tmp/toolchain_sadb_20261119/HISTORY_AFTER.json").write_text(json.dumps(out,ensure_ascii=False,indent=1),encoding="utf8")
E=html.escape
tab=[]
def row(*fields):
    tab.append("<tr>"+"".join("<td>"+s+"</td>" for s in fields)+"</tr>")
for c in cases:
    row(E(c["author"]),E(c["file"])+" · "+str(c["line"]),E(c["before"]),E(c["question"])+"<p>"+E(c["judgment"])+"</p><details><summary>Current source</summary>"+E(c["after"])+"</details>","Applied locally · 3/4")
row("Ben","Chapter 7 · §7.2.2","Isometric Joint Design Factors","<b>Joint Design Factors</b> covers the isometric experiments and the proposed dynamic joint-power measurements.","Applied locally")
row("Ben and print review","Figures 6.2, 6.6, 7.1, 7.3, 7.4, D8","Figure 7.3 arrowheads measured 3.84 × 3.84 pt.","42 enlarged arrowheads across six existing figures. In the final inset layout, Figure 7.3 arrowheads print at 9.44 × 7.55 pt. The Chapter 5 curation figure also has larger arrows.","Applied locally")
row("Independent reader","Figure 5.1 and Figure 5.3","An earlier screenshot and general text made it unclear which collection or control the reader was seeing.","The preserved 943-Paper map is dated relative to the 1,393-entry Table snapshot. The new screenshots show the actual lookup and drill-through. Appendix E explains the denominators.","Applied locally")
row("Ben and PSU print review","Chapters 5–6 · paragraph breaking","A chapter-wide widow and club override and two timing paragraphs locked against breaks were present in the earlier source.","Removed the numeric overrides from the revised chapters. Widow and club penalties control where a paragraph may split across pages. The actual margins, printed folios, font sizes, and line spacing are verified on the built pages.","Applied locally")
row("PSU front-matter review","Title page and copyright","The early title used a different placement grid, and the copyright notice used a month and an additional rights statement.","Checked the title against PSU’s filled dissertation template. The copyright page now reads © 2026 Ben Bolen.","Applied locally with author approval")
row("PSU front-matter review","Abstract and front lists","The Abstract was ii. The contents heading read Contents, and the two lists began below other material and were absent from the printed contents.","The Abstract is i. The Table of Contents, List of Tables, and List of Figures start on their own pages, with 12-point headings. The contents lists the two front lists and their Roman folios.","Applied locally with author approval")
row("Ben and PDF font audit","Chapters 1, 2, 4, 5, and 6, and Appendix C","Twelve active graphics contained 18 unembedded Helvetica or Helvetica-Bold resources.","Embedded the fonts while retaining the original drawing streams, word positions, and all seven SADb figure callouts. The twelve source figures rendered identically at three resolutions in MuPDF.","Applied locally with author approval")
row("Later author edit in the local Overleaf export","Chapter 6 §6.1 · passive toe and suspended motion","The earlier editor capture said ‘I added a passive toe joint.’ The later export changed this to ‘There is a passive toe joint.’","Preserved the later toe-joint wording. Described stepping-like leg movements in air and the fall when support was removed, as Ben requested.","Applied locally")
row("Ben and independent print review","Dissertation · one-inch margins and citation boxes","The original one-point link frames remained visible. Across a page break the PDF also drew extra boxes over an unrelated page header and at the paper edge. Some text extended into the side margins.","Reserved four more PDF points on each side for optical punctuation and citation borders. The page output routine keeps the actual citation boxes, including citations that wrap, and excludes the header from the preceding page’s active link. Both PDFs preserve the visible link borders.","Applied locally at Ben’s direction")
row("Independent final page review","Appendix D · Figure D7 and §D.5.2","At the revised text width, the figure and caption exceeded the available float height by 1.805 TeX points. The shorter caption fit, but the full-page figure split the next subsection’s opening sentence.","Condensed the caption while retaining the 10–12.5 s trial window, support condition, controller identities, line-style definitions, and filtered spiking readouts. Placed Figure D7 before the next subsection. The artwork keeps its legible size.","Applied locally")
row("Independent PSU review","Published Chapter 3 · opening page, source line 5","The coauthored publication is introduced on a page shared with the chapter text. Several contribution assignments differ from the published statement, including data curation.","PSU calls for a separate chapter heading page with a full citation and the published contribution statement. The article credits Bolen and Pang with data curation and credits Bolen and Hunt jointly with several roles.","Open · protected published chapter")
earlier=here.parent/"ch6_advisor_20261119/REVIEW_TABLE.html"
row("Independent figure rendering review","Chapter 4 · Figures 4.9 and 4.12, printed pp. 104 and 113","The incoming PDFs contain two line-dash instructions missing their phase. Poppler stops drawing part of each legend.","Poppler omits the unoptimized and optimized swatches in 4.9 and the six entries from Bend radius through Reference bones in 4.12. Adding the missing zero phase in memory restores both portions. The archive contains the pre-embedding originals.","Open · original Chapter 4 artwork, outside the font and margin changes")
soup=BeautifulSoup(earlier.read_text(encoding="utf8"),"html.parser")
for t in soup.find_all("table"):
    heading=t.find_previous(["summary","h2","h3"])
    heads=[x.get_text(" ",strip=True) for x in t.find_all("th",recursive=True)][:15]
    for tr in t.find_all("tr"):
        if tr.find_parent("table") is not t:continue
        cells=tr.find_all("td",recursive=False)
        if not cells:continue
        vals=[E(c.get_text(" ",strip=True)) for c in cells]
        summary="<br>".join("<b>"+E(heads[k] if k<len(heads) else "Source record")+"</b> "+vals[k] for k in range(2,len(vals)))
        row("Earlier review",E(heading.get_text(" ",strip=True) if heading else "Dissertation review"),vals[0],(vals[1]+"<br>"+summary) if len(vals)>1 else "","As recorded on the review date")
body='<!doctype html><html lang="en"><meta charset="utf-8"><title>Dissertation revision review</title><style>body{font:16px/1.5 system-ui;margin:2rem auto;max-width:1500px;padding:0 1.1rem;color:#1c2c38}table{width:100%;border-collapse:collapse;table-layout:fixed}td,th{padding:.55rem;text-align:left;vertical-align:top;border-bottom:1px solid #cbd4da;overflow-wrap:anywhere}th{background:#e6eef4}td:nth-child(4){width:35%}tr:nth-child(even){background:#f7f9fa}details{margin-top:.5em}a{color:#14557a}</style><h1>Dissertation revision review</h1><p>The first rows cover the latest author corrections and independent findings. The remaining rows retain the earlier advisor comments, tracked changes, author directives, and dispositions. Earlier locations refer to the drafts those readers saw. All 444 attributable advisor observations retain their original aggregate score of −571. The current author directions add no advisor events. A local 3/4 indicates source and contextual verification.</p><table><colgroup><col style="width:13%"><col style="width:14%"><col style="width:23%"><col style="width:39%"><col style="width:11%"></colgroup><thead><tr><th>Source</th><th>Location</th><th>Earlier passage or review item</th><th>Change and rationale</th><th>Disposition</th></tr></thead><tbody>'+''.join(tab)+'</tbody></table></html>'
(here/"REVIEW_TABLE.html").write_text(body,encoding="utf8")
print(len(observed),"actual advisor signals;",out["observed_advisor_score"],"points;",len(cases),"current verified decisions;",len(tab),"combined review rows.")
