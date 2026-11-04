"""Phase 3B MGEdb BLASTN and VF/ARG proximity analysis."""
from pathlib import Path
import csv, json, subprocess, shutil, hashlib
BLAST_FIELDS='qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore'.split()

def build_blastn_command(exe, query, db, out, threads=8, **opts):
    """Build a BLAST+ blastn argv list (outfmt is one complete argument)."""
    fields = "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore"
    cmd = [str(exe), "-query", str(query), "-db", str(db), "-out", str(out), "-outfmt", fields, "-num_threads", str(threads)]
    option_names = {"max_target_seqs": "-max_target_seqs", "evalue": "-evalue", "task": "-task"}
    for key, value in opts.items():
        if value is None or key == "outfmt":
            continue
        flag = option_names.get(key, f"-{key}")
        cmd.extend([flag, str(value)])
    return cmd

def makeblastdb(reference, db_prefix):
 exe=shutil.which('makeblastdb')
 if not exe: raise RuntimeError('makeblastdb executable not found; install blast from bioconda')
 return subprocess.run([exe,'-in',str(reference),'-dbtype','nucl','-out',str(db_prefix)],check=True)

def blastn_search(query,reference,outdir,threads=8,**opts):
 exe=shutil.which('blastn')
 if not exe: raise RuntimeError('blastn executable not found; install blast from bioconda')
 outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); prefix=outdir/'mgedb_blastdb'; makeblastdb(reference,prefix)
 out=outdir/'mge_hits.raw.tsv'; subprocess.run(build_blastn_command(exe,query,prefix,out,threads,**opts),check=True); return out

def _union(intervals):
 merged=[]
 for a,b in sorted((min(a,b),max(a,b)) for a,b in intervals):
  if merged and a<=merged[-1][1]+1: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
  else: merged.append((a,b))
 return merged

def classify(cov, thresholds):
 if cov >= thresholds.get('high_coverage',0.8): return 'HIGH_COVERAGE_REFERENCE_MATCH'
 if cov >= thresholds.get('partial_coverage',0.2): return 'PARTIAL_REFERENCE_MATCH'
 return 'AMBIGUOUS_MATCH'

def parse_hits(raw, metadata, out, thresholds=None, mag_id="UNKNOWN"):
 thresholds=thresholds or {}; md={}
 with open(metadata) as h:
  for r in csv.DictReader(h,delimiter='\t'): md[r.get('record_id','')]=r
 rows=[]
 with open(raw) as h:
  for line in h:
   v=line.rstrip().split('\t')
   if len(v)!=12: continue
   d=dict(zip(BLAST_FIELDS,v)); m=md.get(d['sseqid']); d.update({'mag_id':mag_id, 'contig_id':d['qseqid'], 'query_strand':'+' if int(d['qend'])>=int(d['qstart']) else '-', 'metadata_status':'MAPPED' if m else 'UNMAPPED'})
   if m: d.update({f'mge_{k}':val for k,val in m.items() if k!='record_id'})
   rows.append(d)
 groups={}
 for r in rows: groups.setdefault((r['qseqid'],r['sseqid']),[]).append(r)
 summaries=[]
 for (contig,sid),rs in groups.items():
  merged=_union([(int(r['sstart']),int(r['send'])) for r in rs]); slen=int(rs[0]['slen']); covered=sum(b-a+1 for a,b in merged); cov=covered/max(1,slen)
  summaries.append({'mag_id':rs[0]['mag_id'],'contig_id':contig,'mge_reference_id':sid,'subject_length':slen,'merged_subject_covered_bp':covered,'reference_coverage':f'{cov:.6f}','hsp_count':len(rs),'match_interpretation':classify(cov,thresholds),'metadata_status':rs[0]['metadata_status']})
 fields=sorted({k for r in rows for k in r})
 with open(out,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
 return rows,summaries

def _read(path):
 with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))

def proximity(gene_file, functional_file, mge_file, out, window=5000):
 """Join functional query gene IDs to genes, then compare with MGE query intervals."""
 genes_by_id={r.get('gene_id'):r for r in _read(gene_file) if r.get('gene_id')}
 functional=_read(functional_file); mges=_read(mge_file); seen=set(); outrows=[]
 # One gene is represented once, while each original MGE HSP is retained.
 for fh in functional:
  gid=fh.get('qseqid'); g=genes_by_id.get(gid)
  if not g or gid in seen: continue
  seen.add(gid); ga,gb=sorted((int(g['start']),int(g['end'])))
  for h in mges:
   if g.get('mag_id')!=h.get('mag_id') or g.get('contig_id')!=h.get('contig_id'): continue
   ma,mb=sorted((int(h['qstart']),int(h['qend'])))
   dist=0 if max(ga,ma)<=min(gb,mb) else min(abs(ma-gb),abs(ga-mb))-1
   relation='OVERLAPPING' if dist==0 else ('NEARBY' if dist<=window else 'SAME_CONTIG_DISTANT')
   if relation not in ('OVERLAPPING','NEARBY'): continue
   outrows.append({'mag_id':g['mag_id'],'contig_id':g['contig_id'],'gene_id':gid,'mge_reference_id':h.get('sseqid',''),'gene_start':ga,'gene_end':gb,'mge_query_start':ma,'mge_query_end':mb,'distance_bp':dist,'relation':relation})
 fields=['mag_id','contig_id','gene_id','mge_reference_id','gene_start','gene_end','mge_query_start','mge_query_end','distance_bp','relation']
 with open(out,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(outrows)
 return len(outrows)
