from pathlib import Path
import json

def test_real_qc_counts():
    root=Path(__file__).parents[1]/'build/databases/v1/qc'
    if not root.exists(): return
    assert json.load(open(root/'vfdb_qc.json'))['protein']['records']==4755
    s=json.load(open(root/'sarg_qc.json'))
    assert s['all']['records']==116580
    m=json.load(open(root/'mgedb_qc.json'))
    assert m['nucleotide']['records']==7471 and m['cds']['records']==36338
