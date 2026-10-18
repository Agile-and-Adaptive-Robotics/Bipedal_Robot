# Check PDF text extraction options in the myo env (read-only probe).
libs = ["pdfminer", "PyPDF2", "pypdf", "fitz"]
for l in libs:
    try:
        __import__(l)
        print(l, "OK")
    except Exception as e:
        print(l, "missing")
