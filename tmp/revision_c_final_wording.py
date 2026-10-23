from pathlib import Path
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
p=C/'chapters/15-background.tex';s=p.read_text().replace('I adjusted BPA attachment and via-point locations','I adjusted BPA origins, insertions, and via-point locations');p.write_text(s,encoding='utf-8')
p=C/'chapters/50-database.tex';s=p.read_text().replace('Feedback records describe pathways and link them to supporting papers;','Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers;');p.write_text(s,encoding='utf-8')
p=C/'chapters/30-torque.tex';s=p.read_text().replace('The measured torques meet or exceed the human reference over most of the tested range.','The measured torque magnitudes reach or exceed the human-reference magnitudes over most of the tested range.');p.write_text(s,encoding='utf-8')
print('Made attachment terms, pathway-record example, and torque-magnitude comparison explicit')
