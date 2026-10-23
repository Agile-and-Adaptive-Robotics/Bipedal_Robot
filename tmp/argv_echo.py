import sys
# Echo how the exact /TR argument string arrives after cmd/CRT argument parsing.
for a in sys.argv[1:]:
    print(repr(a))
