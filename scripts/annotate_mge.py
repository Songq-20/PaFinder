from pathlib import Path
import sys,json,shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pafinder.annotation.mge import blastn_search,parse_hits,proximity
import yaml
if len(sys.argv)!=2: raise SystemExit('usage: python scripts/annotate_mge.py MAG.fasta')
root=Path(__file__).resolve().parents[1]; cfg=yaml.safe_load(open(root/'config/phase3b.yaml')); out=root/cfg['output_dir']; out.mkdir(parents=True,exist_ok=True); mag=Path(sys.argv[1]); raw=blastn_search(mag,root/cfg['mge_reference'],out,cfg['threads'],**cfg.get('blastn',{})); rows,sums=parse_hits(raw,root/'build/databases/v1/mgedb/metadata.tsv',out/'mge_hits.tsv',cfg.get('classification',{}),mag.stem)
with open(out/'mge_reference_summary.tsv','w',newline='') as f:
 w=__import__('csv').DictWriter(f,fieldnames=list(sums[0]) if sums else ['mag_id','contig_id','mge_reference_id','reference_coverage','match_interpretation'],delimiter='\t'); w.writeheader(); w.writerows(sums)
for label in ('vf','arg'): proximity(root/cfg['phase3a_dir']/ 'genes.tsv', root/cfg['phase3a_dir']/(label+'_hits.tsv'), out/'mge_hits.tsv', out/(label+'_mge_proximity.tsv'), cfg['neighbor_window'])
json.dump({'mag':str(mag.resolve()),'blast_hsp_count':len(rows),'unique_mge_references':len({r['sseqid'] for r in rows}),'status':'COMPLETED'},open(out/'run_manifest.json','w'),indent=2)
open(out/'mge_summary.tsv','w').write('metric\tvalue\nblast_hsp_count\t'+str(len(rows))+'\nunique_mge_references\t'+str(len({r['sseqid'] for r in rows}))+'\n')
