import sys, urllib.parse, urllib.request, json

def probe(label, value):
    url = ('https://api.openalex.org/works?filter=raw_author_name.search:'
           + urllib.parse.quote(value)
           + ',from_publication_date:2026-09-01,to_publication_date:2026-10-22&per-page=5')
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'SADbMonthlyWatch/1.0'})
        d = json.load(urllib.request.urlopen(req, timeout=30))
        print('%-38s count=%s' % (label, d['meta']['count']))
    except Exception as e:
        print('%-38s ERROR %s' % (label, e))

probe('Frigon', 'Frigon')
probe('Beer|Frigon', 'Beer|Frigon')
probe('Beer|Frigon|Grillner', 'Beer|Frigon|Grillner')
probe('with-Bueschges (ascii)', 'Beer|Buschges|Frigon')
probe('with-Buschges (umlaut)', 'Beer|Büschges|Frigon')
probe('10 ascii surnames', 'Beer|Buschges|Deng|Duysens|Forssberg|Frigon|Grillner|Hunt|Ijspeert|Jankowska')
