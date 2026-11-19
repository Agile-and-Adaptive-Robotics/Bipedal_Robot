"""Apply the approved front-matter and one-inch print settings to the verified baseline."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib,json
from capture_baseline import fingerprint

here=Path(__file__).resolve().parent
source=here.parents[1]/"Final"
pre={"main.tex":(8075,"9c2e1d18"),"psutex.cls":(16896,"bbfa0c9b"),"chapters/01-copyright.tex":(217,"05ba211c")}
backup=here/"FRONT_MATTER_BASELINE.zip"
if not backup.exists():
    candidates={name:(source/name).read_bytes() for name in pre}
    assert all(fingerprint(data.decode("utf8").replace("\r\n","\n"))==pre[name] for name,data in candidates.items())
    with ZipFile(backup,"x",compression=ZIP_DEFLATED) as z:
        for name,data in candidates.items():z.writestr(name,data)
with ZipFile(backup) as z:
    before={name:z.read(name) for name in pre}
strings={name:data.decode("utf8").replace("\r\n","\n") for name,data in before.items()}
main=strings["main.tex"]; cls=strings["psutex.cls"]
def change(s,old,new):
    assert s.count(old)==1,old
    return s.replace(old,new)
main=change(main,r"\usepackage{graphicx,epstopdf}",r"\usepackage{graphicx,epstopdf}"+"\n"+r"\usepackage[pass]{geometry}")
main=change(main,"\n"+r"\usepackage[subfigure]{tocloft}","\n"+r"\usepackage[titles,subfigure]{tocloft}")
main=change(main,r"\maketitle",r"""\hypersetup{pageanchor=false}
\newgeometry{margin=1in,headheight=0pt,headsep=0pt,footskip=0pt}
\maketitle""")
main=change(main,"% Abstract\n",r"""\newgeometry{margin=1in,includehead,headheight=.5in,headsep=.05in}
\pagenumbering{roman}
\hypersetup{pageanchor=true}
\begingroup
\microtypesetup{protrusion=false}
\makeatletter
\renewcommand{\@makeschapterhead}[1]{\noindent{\centering\normalfont\bfseries\normalsize #1\par}\nobreak\vspace{\baselineskip}}
\makeatother

% Abstract
""")
main=change(main,r"\mainmatter",r"""\clearpage
\endgroup
\restoregeometry
\mainmatter""")
main=main.replace(r"\newgeometry{margin=1in,",r"\newgeometry{margin=1in,left=\dimexpr1in+4bp\relax,right=\dimexpr1in+4bp\relax,")
main=change(main,r"\begin{document}",r"""\makeatletter
\ifdefined\pdfrunninglinkoff
  \ifdefined\pdfannot@link@off@@
    \let\pdfannot@link@off@@\pdfrunninglinkoff
    \let\pdfannot@link@on@@\pdfrunninglinkon
    \AddToHook{shipout/before}{\setbox\ShipoutBox=\vbox{\pdfrunninglinkoff\box\ShipoutBox}}
  \fi
\fi
\makeatother
\begin{document}""")
for name in (r"\@lmarwidth",r"\@smarwidth"):
    cls=change(cls,r"\renewcommand{"+name+r"}{0in}",r"\renewcommand{"+name+r"}{4bp}")
cls=change(cls,r"\renewcommand{\@texwidth}{6.5in}",r"\renewcommand{\@texwidth}{\dimexpr6.5in-8bp\relax}")
cls=change(cls,r"  \pagenumbering{roman}",r"  \pagenumbering{gobble}")
start=cls.index(r"\newcommand{\copyrightpage}")
end=cls.index("% Title page",start)
cls=cls[:start]+r"""\newcommand{\copyrightpage}[1][\@commencementyear]{%
  \thispagestyle{empty}%
  \vspace*{2.17in}%
  \begin{center}\normalfont\normalsize\textcopyright~#1\ \@author\end{center}%
  \clearpage
}

"""+cls[end:]
start=cls.index(r"\renewcommand{\titlepage}")
end=cls.index("\n}",start)+2
cls=cls[:start]+r"""\renewcommand{\titlepage}{%
  \thispagestyle{empty}%
  \begingroup
  \centering\normalfont\normalsize
  \setlength{\baselineskip}{13.85pt}%
  \setlength{\parskip}{0pt}%
  \setlength{\topskip}{11.25pt}%
  {\setlength{\baselineskip}{27.7pt}\@title\par}%
  \vspace{5\baselineskip}%
  by\par
  \vspace{\baselineskip}%
  \@author\par
  \vspace{6\baselineskip}%
  A \MakeLowercase{\@doctype} submitted in partial fulfillment of the\\
  requirements for the degree of\par
  \vspace{5\baselineskip}%
  \@graddegree\\in\\\@major\par
  \vspace{5\baselineskip}%
  \@doctype~Committee:\\\@committeemembers\par
  \vspace{6\baselineskip}%
  \@university\\\@commencementyear\par
  \endgroup
  \clearpage
}"""+cls[end:]
strings["main.tex"]=main
strings["psutex.cls"]=cls
strings["chapters/01-copyright.tex"]=r"\copyrightpage[2026]"+"\n"
hashes={}
prior_path=here/"FRONT_MATTER_HASHES.json"
prior=json.loads(prior_path.read_text()) if prior_path.exists() else {}
for name,value in strings.items():
    path=source/name
    local=path.read_bytes()
    eol="\r\n" if b"\r\n" in before[name] else "\n"
    result=value.replace("\n",eol).encode("utf8")
    assert local in (before[name],result) or hashlib.sha256(local).hexdigest()==prior.get(name),name
    if result!=local:path.write_bytes(result)
    hashes[name]=hashlib.sha256(result).hexdigest()
prior_path.write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf8")
