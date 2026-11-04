from pathlib import Path
def test_mgedb_outputs():
 p=Path(__file__).parents[1]/'build/databases/v1/mgedb'
 assert (p/'reference.fna').exists() and (p/'cds_reference.faa').exists()
