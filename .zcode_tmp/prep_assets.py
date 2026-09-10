import zipfile, os, re, shutil, subprocess, sys, json, urllib.request

ROOT = r"D:\GitHub\Bipedal_Robot"
DISS = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation")
TALK = os.path.join(DISS, "Movement and Control of Biomimetic Humanoid Robots.pptx")
OUT  = os.path.join(ROOT, r".zcode_tmp\deck_assets")
os.makedirs(OUT, exist_ok=True)
media_dir = os.path.join(OUT, "talk_media"); os.makedirs(media_dir, exist_ok=True)

# 1) extract talk media + map slide1 (title) image = AARL logo
with zipfile.ZipFile(TALK) as z:
    z.extractall(os.path.join(OUT, "talk_x"))
tx = os.path.join(OUT, "talk_x")
for f in os.listdir(os.path.join(tx, "ppt", "media")):
    shutil.copy(os.path.join(tx, "ppt", "media", f), os.path.join(media_dir, f))
rels = open(os.path.join(tx, "ppt", "slides", "_rels", "slide1.xml.rels")).read()
imgs = re.findall(r'Target="\.\./media/([^"]+)"', rels)
print("TITLE SLIDE IMAGES:", imgs)
for i, im in enumerate(imgs):
    shutil.copy(os.path.join(media_dir, im), os.path.join(OUT, f"aarl_logo_{i}{os.path.splitext(im)[1]}"))

# inventory of media by size for picking
inv = []
for f in sorted(os.listdir(media_dir)):
    p = os.path.join(media_dir, f)
    inv.append((f, os.path.getsize(p)))
print(json.dumps(inv, indent=0))

# 2) get portable ghostscript
gs_dir = os.path.join(OUT, "gs")
gsexe = os.path.join(gs_dir, "gs10.04.0", "bin", "gswin64c.exe")
if not os.path.exists(gsexe):
    req = urllib.request.Request("https://api.github.com/repos/ArtifexSoftware/ghostpdl-downloads/releases",
                                 headers={"User-Agent": "asset-prep"})
    rels_json = json.load(urllib.request.urlopen(req, timeout=60))
    asset_url = None
    for rel in rels_json[:6]:
        for a in rel.get("assets", []):
            if a["name"].endswith("w64.zip"):
                asset_url = a["browser_download_url"]; break
        if asset_url: break
    print("GS asset:", asset_url)
    zp = os.path.join(OUT, "gs.zip")
    urllib.request.urlretrieve(asset_url, zp)
    with zipfile.ZipFile(zp) as z2:
        z2.extractall(gs_dir)
    print("GS extracted:", os.path.exists(gsexe), os.listdir(gs_dir))
print("GSEXE", gsexe)
