import fitz
src = r"D:\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\CAM\05_01_TI_R_008\05_01_TI_R_008_drawing.pdf"
doc = fitz.open(src)
print("pages:", len(doc))
for i, page in enumerate(doc):
    # render at ~300 dpi
    pix = page.get_pixmap(dpi=300)
    out = rf"D:\GitHub\Bipedal_Robot\tmp\sw_cam\drawing_p{i+1}.png"
    pix.save(out)
    print(out, pix.width, "x", pix.height)
    # also dump the text layer
    txt = page.get_text()
    with open(rf"D:\GitHub\Bipedal_Robot\tmp\sw_cam\drawing_text_p{i+1}.txt", "w") as f:
        f.write(txt)
    print("  text chars:", len(txt))
doc.close()
