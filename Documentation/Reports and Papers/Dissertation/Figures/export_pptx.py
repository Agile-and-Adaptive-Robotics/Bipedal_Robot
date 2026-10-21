"""Export the pptx schematics to PDF (dissertation tree) via PowerPoint COM,
and to EPS (components/, for Illustrator) via Ghostscript."""
import os
import win32com.client

FIG = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures"
GS = r"D:\MiKTeX\miktex\bin\x64\mgs.exe"

jobs = [
    (r"10-introduction\components\roadmap.pptx", r"10-introduction\roadmap.pdf"),
    (r"60-sns\components\pipeline.pptx", r"60-sns\pipeline.pdf"),
    (r"60-sns\components\sns_curriculum.pptx", r"60-sns\sns_curriculum.pdf"),
    (r"60-sns\components\sns_library_overview.pptx", r"60-sns\sns_library_overview.pdf"),
    (r"70-overall\components\design_workflow.pptx", r"70-overall\design_workflow.pdf"),
    (r"70-overall\components\sensorium_hierarchy.pptx", r"70-overall\sensorium_hierarchy.pdf"),
    (r"70-overall\components\bpa_construction.pptx", r"70-overall\bpa_construction.pdf"),
]

pp = win32com.client.Dispatch("PowerPoint.Application")
for src, dst in jobs:
    s = os.path.join(FIG, src)
    d = os.path.join(FIG, dst)
    pres = pp.Presentations.Open(s, WithWindow=False)
    pres.SaveAs(d, 32)  # ppSaveAsPDF
    pres.Close()
    eps = os.path.join(FIG, os.path.dirname(src),
                       os.path.basename(dst).replace('.pdf', '.eps'))
    os.system('""%s" -q -dNOPAUSE -dBATCH -sDEVICE=eps2write -sOutputFile="%s" "%s""'
              % (GS, eps, d))
    print("exported", dst)
pp.Quit()
print("done")
