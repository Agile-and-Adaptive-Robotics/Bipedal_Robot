"""Check .zotero-ft-cache availability for the audit's key items."""
import os

STOR = r"C:\Users\Ben Bolen\Zotero\storage"
KEYS = {
    "77QUIHAM": "Renshaw locomotion: antagonism (Hultborn-class)",
    "D3TVUUJP": "IaIN + Renshaw in fictive locomotion",
    "RYXH6JSP": "Renshaw cell activity during fictive locomotion",
    "6MNCI8I6": "Recurrent inhibition of Ia pathway (Hultborn)",
    "F2FNGYXP": "Locomotor-related group Ib pathway (Gossard/Hultborn)",
    "YANK2JE7": "INs mediating reset by extensor group I afferents",
    "X7KBR5U9": "Proprioceptive input resets central rhythm (spinal cat)",
    "MQ4SINGP": "Flexor reflex afferents reset step cycle",
    "9GTNRHKK": "Proprioceptive control of extensor activity (weight support)",
    "M4CHZYBA": "Spinal reflexes Eccles->Lundberg review",
    "R9D8ZJ55": "Spinal control of locomotion cat to man",
}
for k, label in KEYS.items():
    folder = os.path.join(STOR, k)
    if not os.path.isdir(folder):
        print(f"{k}: NO FOLDER")
        continue
    files = os.listdir(folder)
    pdf = [f for f in files if f.lower().endswith(".pdf")]
    cache = ".zotero-ft-cache" in files
    print(f"{k}: {'PDF' if pdf else '-'} {'CACHE' if cache else '-'}  {label}")
