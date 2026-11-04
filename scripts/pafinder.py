import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pafinder.integration.final import run
def main():
 p=argparse.ArgumentParser(prog='pafinder'); sub=p.add_subparsers(dest='cmd'); b=sub.add_parser('bac'); b.add_argument('-i','--input',required=True); b.add_argument('-o','--output',required=True); b.add_argument('--db'); b.add_argument('-t','--threads',type=int,default=8); b.add_argument('--dry-run',action='store_true'); b.add_argument('--plot-colocalization',action='store_true'); a=p.parse_args()
 if a.cmd!='bac': p.error('use: pafinder bac')
 if a.dry_run: print({'phase3a':'build/phase3a','phase3b':'build/phase3b','output':a.output}); raise SystemExit(0)
 run('build/phase3a','build/phase3b',a.output,a.plot_colocalization)
if __name__ == '__main__': main()
