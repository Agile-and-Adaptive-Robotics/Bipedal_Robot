import fitz
from pathlib import Path
r=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\AE_Revision_20261023')
d=fitz.open(r/'output/Bolen_Defense_Assertion_Evidence_20261023_v6.pdf')
print('PDF pages:',len(d))
