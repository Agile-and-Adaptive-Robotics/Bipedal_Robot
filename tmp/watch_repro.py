import sys, urllib.parse
sys.path.insert(0, r'D:\Github\Bipedal_Robot\SADb_audit\watch')
import monthly_watch as m

authors = ['Beer', 'Büschges', 'Deng', 'Duysens', 'Forssberg', 'Frigon', 'Grillner', 'Hunt',
           'Ijspeert', 'Jankowska', 'McCrea', 'Pearson', 'Perreault', 'Prochazka', 'Rossignol',
           'Rybak', 'Shevtsova', 'Szczecinski', 'Young', 'Zhang']
url = ('https://api.openalex.org/works?filter=raw_author_name.search:'
       + urllib.parse.quote('|'.join(authors))
       + ',from_publication_date:2026-09-01,to_publication_date:2026-10-22&per-page=50&cursor=')
for i in range(3):
    try:
        d = m.http_json(url + urllib.parse.quote('*'))
        print('attempt', i, 'OK count', d['meta']['count'])
        break
    except RuntimeError as e:
        print('attempt', i, 'FAIL', str(e)[:200])
