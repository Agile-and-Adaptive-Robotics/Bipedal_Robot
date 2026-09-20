"""Check: is recXVuB69yaIn6437 (Duysens & Forner-Cordero 2018) already linked to a Review record?"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
records = json.load(open(os.path.join(os.path.dirname(HERE), "export", "sadb_export.json"), encoding="utf-8"))
for r in records:
    if r["id"] == "recXVuB69yaIn6437":
        print("reviews:", r["reviews"], "| has_notes:", r["has_notes"])
# any record whose reviews mention Duysens?
for r in records:
    if any("Duysens" in x for x in r["reviews"]):
        print("linked:", r["id"], "|", r["title"][:60], "->", r["reviews"])
