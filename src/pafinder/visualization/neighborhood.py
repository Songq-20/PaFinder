from pathlib import Path
import csv,json

def read_tsv(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def _short(x):
 x=x or ''; return x if len(x)<=24 else x[:21]+'...'
def _arrow(ax,start,end,y,strand,color,label,edge='#333333',alpha=1):
 import matplotlib.patches as patches
 h=.28; head=min(max((end-start)*.18,20),max(40,(end-start)))
 if strand=='-': pts=[(end,y-h),(start+head,y-h),(start,y),(start+head,y+h),(end,y+h)]
 elif strand=='+': pts=[(start,y-h),(end-head,y-h),(end,y),(end-head,y+h),(start,y+h)]
 else: pts=[(start,y-h),(end,y-h),(end,y+h),(start,y+h)]
 ax.add_patch(patches.Polygon(pts,closed=True,facecolor=color,edgecolor=edge,linewidth=1,alpha=alpha))
 if label: ax.text((start+end)/2,y+h+.04,label,ha='center',va='bottom',fontsize=6,rotation=45)

def plot_regions(phase3a,phase3b,outdir,flank=5000,max_plots=10,formats=('svg','png')):
 import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
 out=Path(outdir); out.mkdir(parents=True,exist_ok=True); genes=read_tsv(Path(phase3a)/'genes.tsv'); vfh=read_tsv(Path(phase3a)/'vf_hits.tsv'); argh=read_tsv(Path(phase3a)/'arg_hits.tsv'); mges=read_tsv(Path(phase3b)/'mge_hits.tsv'); prox=read_tsv(Path(phase3b)/'vf_mge_proximity.tsv')+read_tsv(Path(phase3b)/'arg_mge_proximity.tsv')
 gby={g['gene_id']:g for g in genes}; types={}
 for r in vfh: types.setdefault(r['qseqid'],set()).add('VF')
 for r in argh: types.setdefault(r['qseqid'],set()).add('ARG')
 # unique genomic anchors, prioritize overlap then distance
 anchors={}
 for r in prox:
  if r['relation'] not in ('OVERLAPPING','NEARBY') or r['gene_id'] not in gby: continue
  key=(r['mag_id'],r['contig_id'],gby[r['gene_id']]['start'],gby[r['gene_id']]['end'])
  anchors.setdefault(key,[]).append(r)
 ordered=sorted(anchors.items(),key=lambda kv:(0 if any(x['relation']=='OVERLAPPING' for x in kv[1]) else 1,min(int(x['distance_bp']) for x in kv[1])))
 if max_plots!='all': ordered=ordered[:int(max_plots)]
 index=[]; qc={'candidate_regions':len(anchors),'plots':0,'missing_genes':0}
 for i,((mag,contig,gs,ge),prs) in enumerate(ordered,1):
  center=(int(gs)+int(ge))//2; lo=max(1,center-flank); hi=center+flank; locus=[g for g in genes if g['mag_id']==mag and g['contig_id']==contig and int(g['end'])>=lo and int(g['start'])<=hi]
  mh=[m for m in mges if m['mag_id']==mag and m['contig_id']==contig and int(m['qend'])>=lo and int(m['qstart'])<=hi]
  fig,ax=plt.subplots(figsize=(6.3,2.6)); ax.set_xlim(lo,hi); ax.set_ylim(-.8,1.3); ax.set_yticks([0,1]); ax.set_yticklabels(['ORFs','MGE reference alignment']); ax.set_title('VF/ARG–MGE homologous genomic neighborhood',fontsize=10)
  for m in mh:
   a,b=sorted((int(m['qstart']),int(m['qend']))); label=_short(m.get('mge_Element_name') or m.get('mge_category') or 'MGE reference alignment'); ax.add_patch(__import__('matplotlib').patches.Rectangle((a,.82),b-a,.25,facecolor='#D99A56',edgecolor='#8a5a20',linestyle='--',alpha=.45)); ax.text((a+b)/2,1.1,label,ha='center',fontsize=6)
  for g in locus:
   a,b=sorted((int(g['start']),int(g['end']))); ts=types.get(g['gene_id'],set()); color='#D87C76' if ts=={'VF'} else '#E1AC48' if ts=={'ARG'} else '#B8BEC7'; edge='#6b3b8f' if ts=={'VF','ARG'} else '#333333'; _arrow(ax,a,b,0,g.get('strand',''),color,g['gene_id'].split('__')[-1],edge)
  ax.set_xlabel(f'Genomic position on {contig} (1-based)'); ax.grid(axis='x',alpha=.2); fig.tight_layout()
  stem=f'locus_{i:03d}'
  for fmt in formats: fig.savefig(out/(stem+'.'+fmt),dpi=300 if fmt=='png' else None,facecolor='white')
  plt.close(fig); index.append({'plot_id':stem,'mag_id':mag,'contig_id':contig,'window_start':lo,'window_end':hi,'vf_arg_genes':','.join(sorted({x['gene_id'] for x in prs})),'mge_references':','.join(sorted({x['mge_reference_id'] for x in prs})),'proximity_rows':len(prs)})
 qc['plots']=len(index)
 with open(out/'plot_index.tsv','w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(index[0]) if index else ['plot_id','mag_id','contig_id'],delimiter='\t'); w.writeheader(); w.writerows(index)
 json.dump(qc,open(out/'plot_qc.json','w'),indent=2); return qc
