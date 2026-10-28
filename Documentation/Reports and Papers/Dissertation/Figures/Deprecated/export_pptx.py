"""Export editable dissertation schematics to PDF and EPS using PowerPoint.

Run make_pptx_figures.py first. The manifest limits exports to authored
schematics; experimental plots and the isometric-force chapter are excluded.
"""
from pathlib import Path
import argparse, json, subprocess
import win32com.client

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--figure-root',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();root=args.figure_root
    jobs=json.loads((root/'components/schematic_jobs.json').read_text())
    gs=Path(r'D:\MiKTeX\miktex\bin\x64\mgs.exe')
    app=win32com.client.DispatchEx('PowerPoint.Application')
    try:
        for chapter,name in jobs:
            src=root/chapter/'components'/(name+'.pptx')
            dst=root/chapter/(name+'.pdf')
            doc=app.Presentations.Open(str(src),ReadOnly=True,WithWindow=False)
            try:doc.SaveAs(str(dst),32)
            finally:doc.Close()
            subprocess.run([str(gs),'-q','-dNOPAUSE','-dBATCH','-sDEVICE=eps2write',f'-sOutputFile={src.with_suffix(".eps")}',str(dst)],check=True)
            print('Exported',chapter,name,flush=True)
    finally:app.Quit()

if __name__=='__main__':main()
