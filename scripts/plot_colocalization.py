import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pafinder.visualization.neighborhood import plot_regions
p=argparse.ArgumentParser(); p.add_argument('--phase3a',required=True); p.add_argument('--phase3b',required=True); p.add_argument('-o','--output',required=True); p.add_argument('--flank',type=int,default=5000); p.add_argument('--max-plots',default='10'); p.add_argument('--format',default='svg,png'); a=p.parse_args(); print(plot_regions(a.phase3a,a.phase3b,Path(a.output)/'colocalization',a.flank,a.max_plots,tuple(a.format.split(','))))
