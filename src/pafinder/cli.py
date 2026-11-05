import argparse
from pathlib import Path
from pafinder.pipeline.bacteria import run_pipeline

def main():
 p=argparse.ArgumentParser(prog='pafinder'); s=p.add_subparsers(dest='cmd'); b=s.add_parser('bac'); b.add_argument('-i','--input',required=True); b.add_argument('-o','--output',required=True); b.add_argument('--db',required=True); b.add_argument('-t','--threads',type=int,default=8); b.add_argument('--dry-run',action='store_true'); b.add_argument('--plot-colocalization',action='store_true'); b.add_argument('--skip-skani',action='store_true'); a=p.parse_args()
 if a.cmd!='bac': p.error('use pafinder bac')
 statuses=run_pipeline(a.input,a.output,a.db,a.threads,a.plot_colocalization,dry_run=a.dry_run,skip_skani=a.skip_skani)
 if any(x.get('status')=='FAILED' for x in statuses): raise SystemExit(1)
if __name__=='__main__': main()
