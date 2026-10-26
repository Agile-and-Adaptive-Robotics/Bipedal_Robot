# Read-only peek at ZCode's tasks-index.sqlite: table names + row counts only.
import sqlite3

p = r"C:\Users\Ben Bolen\.zcode\v2\tasks-index.sqlite"
con = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
cur = con.cursor()
tables = [r[0] for r in cur.execute(
    "SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
for t in tables:
    n = cur.execute(f'SELECT COUNT(*) FROM "{t}"').fetchone()[0]
    cols = [c[1] for c in cur.execute(f'PRAGMA table_info("{t}")').fetchall()]
    print(t, "| rows:", n, "| cols:", cols[:12])
con.close()
