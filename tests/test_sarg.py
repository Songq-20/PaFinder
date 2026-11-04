from pathlib import Path
def test_sarg_sources():
 p=Path(__file__).parents[1]/'build/databases/v1/qc/sarg_qc.json'
 assert p.exists()
