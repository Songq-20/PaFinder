from pathlib import Path
from pafinder.pipeline.bacteria import run_pipeline

def test_pipeline_dry_run(tmp_path):
 inp=tmp_path/'mags'; inp.mkdir(); (inp/'a.fna').write_text('>c\nATG\n'); out=tmp_path/'out'
 r=run_pipeline(inp,out,Path('.'),dry_run=True)
 assert r[0]['status']=='NOT_RUN'; assert (out/'qc'/'mag_status.tsv').exists()
