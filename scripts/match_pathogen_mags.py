import argparse,csv,json,subprocess,shutil
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--db',required=True); p.add_argument('--manifest',required=True); p.add_argument('--output',required=True); p.add_argument('--threads',type=int,default=8); a=p.parse_args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
root=Path(a.input); mags=sorted([root] if root.is_file() else [x for x in root.iterdir() if x.suffix.lower() in ('.fa','.fna','.fasta','.gz')]); raw=out/'skani_hits.raw.tsv'; exe=shutil.which('skani')
if not exe: raise SystemExit('skani not found')
subprocess.run([exe,'search','-d',a.db,'-t',str(a.threads),'-o',str(raw),*map(str,mags)],check=True)
json.dump({'query_count':len(mags),'threads':a.threads,'status':'COMPLETED'},open(out/'run_manifest.json','w'),indent=2)
