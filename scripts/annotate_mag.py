from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pafinder.annotation.mag import run
if len(sys.argv)!=2: raise SystemExit('usage: python scripts/annotate_mag.py MAG.fasta')
run(Path(__file__).resolve().parents[1]/'config/phase3a.yaml',sys.argv[1])
