import csv

rows = list(csv.DictReader(open(r'SADb_audit\export\sadb_export.csv', encoding='utf-8')))
live = [r for r in rows if r['archived'].strip().lower() not in ('true', '1')]
arch = [r for r in rows if r['archived'].strip().lower() in ('true', '1')]

def nz(v):
    return v is not None and v.strip() != ''

print('total rows', len(rows), 'live', len(live), 'archived', len(arch))
for name, col in [('notes', 'notes'), ('has_notes_flag', 'has_notes'),
                  ('afferents', 'afferents'), ('robot_sim', 'robot_sim'),
                  ('feedback', 'feedback'), ('animals', 'animals'),
                  ('has_pdf', 'has_pdf'), ('doi', 'doi')]:
    print(name, 'live-filled', sum(1 for r in live if nz(r[col])),
          'all-filled', sum(1 for r in rows if nz(r[col])))
