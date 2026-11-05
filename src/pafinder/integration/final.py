from pathlib import Path
import csv,json

def read(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(rows,p,fields=None):
 p.parent.mkdir(parents=True,exist_ok=True); fields=fields or (list(rows[0]) if rows else ['status'])
 with open(p,'w',newline='') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)

def cluster_mge_loci(mge_rows, proximity_rows):
 groups={}
 for r in mge_rows:
  a,b=sorted((int(r['qstart']),int(r['qend']))); groups.setdefault((r['mag_id'],r['contig_id']),[]).append((a,b,r))
 loci=[]
 for (mag,contig),items in groups.items():
  items.sort(key=lambda x:(x[0],x[1])); cur=None
  for a,b,r in items:
   if cur is None or a>cur['end']+0:
    if cur:loci.append(cur)
    cur={'mag_id':mag,'contig_id':contig,'start':a,'end':b,'references':{r.get('sseqid','')},'hsp_count':1}
   else:
    cur['end']=max(cur['end'],b); cur['references'].add(r.get('sseqid','')); cur['hsp_count']+=1
  if cur:loci.append(cur)
 for i,l in enumerate(loci,1):
  genes={r['gene_id'] for r in proximity_rows if r['mag_id']==l['mag_id'] and r['contig_id']==l['contig_id'] and int(r['mge_query_end'])>=l['start'] and int(r['mge_query_start'])<=l['end']}; l.update({'locus_id':f"{l['mag_id']}__{l['contig_id']}__mge_locus_{i:04d}",'reference_ids':';'.join(sorted(l.pop('references'))),'associated_gene_ids':';'.join(sorted(genes))})
 return loci

def run(phase3a,phase3b,out,plot=False):
 a,b,o=Path(phase3a),Path(phase3b),Path(out); ev=o/'evidence'; qc=o/'qc'; ev.mkdir(parents=True,exist_ok=True); qc.mkdir(exist_ok=True)
 genes=read(a/'genes.tsv'); vf=read(a/'vf_hits.tsv'); arg=read(a/'arg_hits.tsv'); mge=read(b/'mge_hits.tsv'); vp=read(b/'vf_mge_proximity.tsv'); ap=read(b/'arg_mge_proximity.tsv')
 for src,dst in [(vf,ev/'vf_hits.tsv'),(arg,ev/'arg_hits.tsv'),(mge,ev/'mge_hits.tsv'),(vp,ev/'vf_mge_proximity.tsv'),(ap,ev/'arg_mge_proximity.tsv')]: write(src,dst)
 loci=cluster_mge_loci(mge,vp+ap); write(loci,ev/'mge_loci.tsv',['locus_id','mag_id','contig_id','start','end','reference_ids','hsp_count','associated_gene_ids'])
 vfgenes={r['qseqid'] for r in vf}; arggenes={r['qseqid'] for r in arg}; mags=sorted({g['mag_id'] for g in genes})
 gene_summary=[]
 for g in genes:
  gid=g['gene_id']; flags=[]
  if gid in vfgenes: flags.append('REVIEW_VF_FUNCTION')
  if gid in arggenes: flags.append('REVIEW_ARG_FUNCTION')
  if any(r['gene_id']==gid for r in vp+ap): flags.append('INSPECT_MGE_CONTEXT')
  gene_summary.append({**g,'vf_homolog':str(gid in vfgenes),'arg_homolog':str(gid in arggenes),'evidence_flags':';'.join(flags),'interpretation_status':'HOMOLOGY_REQUIRES_REVIEW' if flags else 'NO_RELEVANT_HOMOLOGY'})
 write(gene_summary,ev/'gene_summary.tsv')
 mags_rows=[]
 for mag in mags:
  gs=[g for g in genes if g['mag_id']==mag]; gids={g['gene_id'] for g in gs}; v=[r for r in vf if r['qseqid'] in gids]; ar=[r for r in arg if r['qseqid'] in gids]; mh=[r for r in mge if r['mag_id']==mag]; vpr=[r for r in vp if r['mag_id']==mag]; apr=[r for r in ap if r['mag_id']==mag]; flags=[]
  if v: flags.append('REVIEW_VF_FUNCTION')
  if ar: flags.append('REVIEW_ARG_FUNCTION')
  if vpr or apr: flags.append('INSPECT_MGE_CONTEXT')
  mags_rows.append({'mag_id':mag,'predicted_gene_count':len(gs),'vf_raw_homolog_genes':len({r['qseqid'] for r in v}),'vf_raw_alignment_rows':len(v),'arg_raw_homolog_genes':len({r['qseqid'] for r in ar}),'arg_raw_alignment_rows':len(ar),'mge_hsp_count':len(mh),'mge_unique_reference_count':len({r['sseqid'] for r in mh}),'vf_mge_pair_count':len(vpr),'arg_mge_pair_count':len(apr),'vf_mge_overlapping':sum(r['relation']=='OVERLAPPING' for r in vpr),'vf_mge_nearby':sum(r['relation']=='NEARBY' for r in vpr),'arg_mge_overlapping':sum(r['relation']=='OVERLAPPING' for r in apr),'arg_mge_nearby':sum(r['relation']=='NEARBY' for r in apr),'evidence_flags':';'.join(flags),'interpretation_status':'SUCCESS_WITH_HITS' if (v or ar or mh) else 'SUCCESS_NO_HITS'})
 write(mags_rows,o/'summary.tsv'); write([{'metric':'mag_count','value':len(mags)},{'metric':'gene_count','value':len(genes)},{'metric':'vf_raw_homolog_genes','value':len(vfgenes)},{'metric':'arg_raw_homolog_genes','value':len(arggenes)},{'metric':'mge_hsp_count','value':len(mge)},{'metric':'vf_mge_pairs','value':len(vp)},{'metric':'arg_mge_pairs','value':len(ap)}],o/'summary_stats.tsv',['metric','value']); write([{'mag_id':m,'status':'SUCCESS_WITH_HITS'} for m in mags],qc/'mag_status.tsv',['mag_id','status']); json.dump({'status':'SUCCESS_WITH_HITS','mag_count':len(mags),'gene_count':len(genes),'vf_alignment_rows':len(vf),'arg_alignment_rows':len(arg),'mge_hsp_count':len(mge)},open(qc/'run_manifest.json','w'),indent=2); json.dump({'phase3a':'existing outputs','phase3b':'existing outputs'},open(qc/'database_versions.json','w'),indent=2)
 if plot:
  try:
   from pafinder.visualization.neighborhood import plot_regions
   plot_regions(a,b,o/'figures',5000,10,('svg','png'))
  except ImportError as e:
   (o/'plot_error.txt').write_text('matplotlib required for plotting: '+str(e))
 return mags_rows
