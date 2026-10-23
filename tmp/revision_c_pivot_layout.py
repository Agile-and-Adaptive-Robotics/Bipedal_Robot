from pathlib import Path
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
p=C/'chapters/50-database.tex';s=p.read_text()
s=s.replace(r'\includegraphics[width=0.96\textwidth]{sadb_pivot_app.png}',r'\includegraphics[width=0.96\textwidth,trim=0 610bp 980bp 150bp,clip]{sadb_pivot_app.png}')
s=s.replace('Application pivot view, with afferent classifications as rows and animal preparations as columns.','Excerpt of the application pivot matrix, with afferent classifications as rows and selected animal preparations as columns.')
p.write_text(s,encoding='utf-8')
p=C/'review/03_SADb.md';s=p.read_text();s+='The pivot screenshot is cropped through LaTeX to show the matrix excerpt at a readable scale and remove the empty application viewport; the underlying PNG is unchanged.\n';p.write_text(s,encoding='utf-8')
print('Removed empty pivot viewport through LaTeX figure framing')
