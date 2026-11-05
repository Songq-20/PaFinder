from pathlib import Path
from pafinder.pipeline.bacteria import run_pipeline

def test_pipeline_dry_run(tmp_path):
 inp=tmp_path/'mags'; inp.mkdir(); (inp/'a.fna').write_text('>c\nATG\n'); out=tmp_path/'out'
 r=run_pipeline(inp,out,Path('.'),dry_run=True)
 assert r[0]['status']=='NOT_RUN'; assert (out/'qc'/'mag_status.tsv').exists()

def test_two_mag_output_layout_isolated(tmp_path):
 inp=tmp_path/'mags'; inp.mkdir()
 for name in ('MAG_001','MAG_002'):
  (inp/(name+'.fna')).write_text('>c\nATG\n')
 out=tmp_path/'results'; from pafinder.pipeline.bacteria import run_pipeline
 statuses=run_pipeline(inp,out,tmp_path,dry_run=True)
 assert len(statuses)==2
 assert (out/'qc'/'mag_status.tsv').exists()
 assert {p.name for p in (out/'mags').iterdir()} == {'MAG_001','MAG_002'}
 assert not (out/'MAG_001').exists() and not (out/'MAG_002').exists()
 assert set(p.name for p in out.iterdir()) <= {'mags','qc'}
