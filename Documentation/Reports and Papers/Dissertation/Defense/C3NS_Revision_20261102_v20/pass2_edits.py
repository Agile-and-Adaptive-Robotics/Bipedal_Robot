from pathlib import Path
W=Path(__file__).resolve().parent
p=W/'revise.mjs';s=p.read_text(encoding='utf-8')
def rep(a,b):
 global s
 assert a in s,a
 s=s.replace(a,b)
rep("if(crop&&ext==='.png'){", "if(crop&&(ext==='.png'||ext==='.gif')){")
rep("const iw=blob.readUInt32BE(16)*(1-crop.left-crop.right),ih=blob.readUInt32BE(20)*(1-crop.top-crop.bottom)", "const iw=(ext==='.png'?blob.readUInt32BE(16):blob.readUInt16LE(6))*(1-crop.left-crop.right),ih=(ext==='.png'?blob.readUInt32BE(20):blob.readUInt16LE(8))*(1-crop.top-crop.bottom)")
rep("tx(5,'Sensory feedback',112,345,250,36,23,false,INK);","tx(5,'Sensory\\nfeedback',35,337,175,62,23,false,INK,'center');")
rep("im(n,a('knee_cutout_approved.png'),455,120,445,230)","im(n,a('knee_cutout_approved.png'),455,169,445,191)")
rep("im(n,old(file),170,104,620,300);im(n,a(eq),n===20?185:100,423,n===20?590:760,46);","im(n,old(file),150,104,660,320);im(n,a(eq),n===20?200:135,435,n===20?560:690,34);")
rep('The force path determines the moment arm about the joint','Muscle geometry determines the joint’s moment arm')
rep("'l₀; lten (cm)'","'l₀; lₜₑₙ (cm)'")
rep("230,118,310,35","229,118,348,35")
rep("570,118,325,35","577,118,348,35")
start=s.index('// Substitute the requested actual pinned-extensor photograph;')
end=s.index("begin(36,",start)
s=s[:start]+s[end:]
rep("const order=[", "await (await import('./pass2_content.mjs')).apply({d,rec,begin,newSlide,im,tx,line,block,a,old,fig,changed,GREEN,INK,W});\nconst order=[")
rep('124,125,38,39','124,125,129,38,39')
p.write_text(s,encoding='utf-8')
p=W/'package_deck.py';s=p.read_text(encoding='utf-8')
s=s.replace('len(order)==128','len(order)==129').replace('range(1,129)','range(1,130)').replace('range(124,129)','range(124,130)').replace("'main_slides':75","'main_slides':76").replace("'total_slides':128","'total_slides':129")
p.write_text(s,encoding='utf-8')
print('Pass 2 authoring and package counts updated')
