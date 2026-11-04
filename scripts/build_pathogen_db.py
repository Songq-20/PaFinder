import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pafinder.databases.pathogen import build_manifest,validate_references
p=argparse.ArgumentParser(); p.add_argument('--metadata'); p.add_argument('--file-list'); p.add_argument('--audit'); p.add_argument('--metadata-only',action='store_true'); p.add_argument('--manifest'); p.add_argument('--genomes-dir'); p.add_argument('--build-skani',action='store_true'); p.add_argument('--threads',type=int,default=8); p.add_argument('--output-dir',default='build/databases/v1/pathogen_references'); a=p.parse_args()
if a.metadata_only or a.manifest is None:
 build_manifest(a.metadata or 'PaFinderDB/pathoMAG_metadata.csv',a.file_list or 'PaFinderDB/pathoMAGid.txt',a.audit,a.output_dir)
else:
 if not a.genomes_dir: p.error('--genomes-dir is required for validation')
 print(validate_references(a.manifest,a.genomes_dir,a.output_dir,a.threads,a.build_skani))
