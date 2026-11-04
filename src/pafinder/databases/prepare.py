from pathlib import Path
import csv,json
try:
 import yaml
except ImportError:
 yaml=None
from .vfdb import prepare as vf
from .sarg import prepare as sg
from .mgedb import prepare as mg
from .common import sha256

def run(config):
 root=Path(__file__).resolve().parents[3]; cfg=yaml.safe_load(open(config)) if yaml else {'input_root':'./PaFinderDB','output_root':'./build/databases/v1'}; inp=(root/cfg['input_root']).resolve(); out=(root/cfg['output_root']).resolve(); out.mkdir(parents=True,exist_ok=True)
 vr,vqc,vmap,vev=vf(inp,out,inp/'PaFinder_VFDB_v1'); sr,sqc,smap=sg(inp,out); mr,mqc,mmap=mg(inp,out)
 def tab(rows,path):
  path.parent.mkdir(parents=True,exist_ok=True)
  with open(path,'w',newline='') as h:
   w=csv.DictWriter(h,fieldnames=sorted({k for r in rows for k in r}),delimiter='\t'); w.writeheader(); w.writerows(rows)
 tab(vr,out/'vfdb/metadata.tsv'); tab(vmap,out/'vfdb/id_mapping.tsv'); tab(vev,out/'vfdb/evidence_review.tsv');
 # Family-level researcher summary is retained as an auditable auxiliary table.
 fam=inp/'PaFinder_VFDB_v1/vfdb_family_summary.tsv'
 if fam.exists():
  import shutil; shutil.copyfile(fam,out/'vfdb/family_summary.tsv')
 tab(sr,out/'sarg_plus/metadata.tsv'); tab(smap,out/'sarg_plus/id_mapping.tsv'); tab(mr,out/'mgedb/metadata.tsv'); tab(mmap,out/'mgedb/id_mapping.tsv')
 qc=out/'qc'; qc.mkdir(exist_ok=True)
 for n,d in [('vfdb',vqc),('sarg',sqc),('mgedb',mqc)]: json.dump(d,open(qc/f'{n}_qc.json','w'),indent=2)
 manifest={'program_version':'0.1.0','input_root':str(inp),'databases':{str(p.relative_to(root)):sha256(p) for p in inp.rglob('*') if p.is_file()}}
 json.dump(manifest,open(out/'manifest.json','w'),indent=2)
 open(qc/'database_summary.md','w').write('# PaFinder database summary\n\n'+json.dumps({'vfdb':vqc,'sarg':sqc,'mgedb':mqc},indent=2))
