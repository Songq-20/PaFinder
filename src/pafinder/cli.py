import argparse
from pathlib import Path
from pafinder.integration.final import run
def main():
 p=argparse.ArgumentParser(prog='pafinder'); sub=p.add_subparsers(dest='cmd'); b=sub.add_parser('bac'); b.add_argument('-i','--input',required=True); b.add_argument('-o','--output',required=True); b.add_argument('--db'); b.add_argument('-t','--threads',type=int,default=8); b.add_argument('--dry-run',action='store_true'); b.add_argument('--plot-colocalization',action='store_true'); a=p.parse_args()
 if a.cmd!='bac': p.error('use: pafinder bac')
 if a.dry_run: print({'phase3a':str(Path(a.db or '.').resolve()/'build/phase3a'),'phase3b':str(Path(a.db or '.').resolve()/'build/phase3b'),'output':a.output}); raise SystemExit(0)
 phase3a=Path(a.db or '.').resolve()/'build/phase3a'
 phase3b=Path(a.db or '.').resolve()/'build/phase3b'
 if not phase3a.exists() or not phase3b.exists():
  raise SystemExit('Existing Phase 3A/3B outputs are required; full MAG execution is not yet wired into this release')
 run(phase3a,phase3b,a.output,a.plot_colocalization)
if __name__ == '__main__': main()
