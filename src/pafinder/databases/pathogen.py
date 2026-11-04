"""Phase 3C pathogen reference audit, skani indexing and matching helpers."""
from pathlib import Path
import csv,json,hashlib,subprocess,shutil,time

def _read_meta(path):
 with open(path,encoding='utf-8-sig',newline='') as h:
  rows=list(csv.reader(h)); header=rows[1] if rows and rows[0] and rows[0][0] != 'Genus_species' else rows[0]; data=rows[2:] if header is rows[1] else rows[1:]
  return [dict(zip(header,r)) for r in data if any(r)]
def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def fasta_valid(p):
 try:
  n=0; has=False
  for line in open(p,errors='replace'):
   if line.startswith('>'): n+=1; has=True
  return has and n>0
 except OSError:return False

def build_manifest(metadata,file_list=None,audit=None,outdir='build/databases/v1/pathogen_references'):
 out=Path(outdir); out.mkdir(parents=True,exist_ok=True); meta=_read_meta(metadata); names=[]
 audit_rows=[]
 if audit and Path(audit).exists():
  ap=Path(audit)/'pathoMAG_full_audit.tsv'
  if ap.exists(): audit_rows=list(csv.DictReader(open(ap),delimiter='\t'))
 audit_by={r.get('assembly_accession',''):r for r in audit_rows}
 if file_list and Path(file_list).exists(): names=[Path(x.strip()).name for x in open(file_list) if x.strip()]
 rows=[]; seen=set()
 for i,r in enumerate(meta):
  acc=str(r.get('accession','')).strip() or str(r.get('match id','')).strip() or f'UNIDENTIFIED_{i+1:04d}'
  if acc in seen: continue
  seen.add(acc); stem=str(r.get('File name','')).strip(); fname=next((n for n in names if acc in n or (stem and stem in n)), stem+'.fna' if stem else acc+'.fna')
  ar=audit_by.get(acc,{}); priority=ar.get('audit_priority','AUDIT_NOT_FOUND'); status='ACCEPTED' if priority in ('NO_TAXONOMIC_LABEL_CONFLICT_DETECTED','GTDB_SUFFIX_CONVENTION') else ('EXCLUDED' if priority in ('CRITICAL_LABEL_CONFLICT','CRITICAL_SUBSTITUTE_REFERENCE') else 'REVIEW_REQUIRED')
  rows.append({'reference_id':acc,'assembly_accession':acc,'fasta_filename':fname,'fasta_path':'NOT_CHECKED','original_species':r.get('Genus_species','UNKNOWN'),'ncbi_species':r.get('ncbi_organism_name','UNKNOWN'),'gtdb_species':r.get('Species','UNKNOWN'),'completeness':r.get('Completeness','UNKNOWN'),'contamination':r.get('Contamination','UNKNOWN'),'reference_status':status,'sequence_check_status':'NOT_CHECKED','pathogenicity_evidence_scope':'metadata candidate; not independently validated','evidence_source':'pathoMAG_metadata.csv','review_notes':ar.get('audit_notes','') or ('audit not found' if not ar else priority),'audit_priority':priority,'taxon_relation':ar.get('taxon_relation','UNKNOWN'),'metadata_review_status':status})
 fields=list(rows[0]) if rows else ['reference_id','assembly_accession','fasta_filename','reference_status']
 for name,subset in [('reference_manifest.tsv',rows),('references_accepted.tsv',[r for r in rows if r['reference_status']=='ACCEPTED']),('references_review_required.tsv',[r for r in rows if r['reference_status']=='REVIEW_REQUIRED']),('references_excluded.tsv',[r for r in rows if r['reference_status']=='EXCLUDED'])]:
  with open(out/name,'w',newline='') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(subset)
 json.dump({'metadata_records':len(rows),'accepted':sum(r['reference_status']=='ACCEPTED' for r in rows),'review_required':sum(r['reference_status']=='REVIEW_REQUIRED' for r in rows),'excluded':sum(r['reference_status']=='EXCLUDED' for r in rows),'sequence_checks':'NOT_CHECKED'},open(out/'database_qc_report.json','w'),indent=2)
 return rows

def build_skani_command(reference_list,out_prefix,threads=8,exe='skani'): return [str(exe),'sketch','-l',str(reference_list),'-o',str(out_prefix),'-t',str(threads)]
def search_skani_command(db,queries,out,threads=8): return ['skani','search','-d',str(db),'-t',str(threads),'-o',str(out),*map(str,queries)]

def parse_skani(path,manifest,thresholds=None):
 thresholds=thresholds or {'ani':95,'af_reference':50,'af_query':50}; md={r['reference_id']:r for r in csv.DictReader(open(manifest),delimiter='\t')}; rows=[]
 with open(path) as h:
  for r in csv.DictReader(h,delimiter='\t'):
   ref=r.get('Ref_file',''); key=Path(ref).stem; m=md.get(key) or next((x for x in md.values() if x.get('fasta_path')==ref),None); d={'query_mag_id':Path(r.get('Query_file','')).stem,'reference_id':m['reference_id'] if m else key,'ANI':r.get('ANI',''),'AF_reference':r.get('Align_fraction_ref',''),'AF_query':r.get('Align_fraction_query',''),'metadata_status':'MAPPED' if m else 'UNMAPPED'}
   if m:d.update(m)
   try:d['match_status']='SPECIES_LEVEL_CANDIDATE' if float(d['ANI'])>=thresholds['ani'] and float(d['AF_reference'])>=thresholds['af_reference'] and float(d['AF_query'])>=thresholds['af_query'] else 'RELATED_REFERENCE_BELOW_THRESHOLD'
   except ValueError:d['match_status']='AMBIGUOUS_SPECIES_MATCH'
   rows.append(d)
 return rows

def validate_references(manifest, genomes_dir, outdir=None, threads=8, build_skani=False):
 """Linux-side validation; resolves portable fasta_filename against genomes_dir."""
 manifest=Path(manifest); out=Path(outdir or manifest.parent); root=Path(genomes_dir).resolve(); rows=list(csv.DictReader(open(manifest),delimiter='\t')); valid=[]; issues=[]
 for r in rows:
  p=root/r['fasta_filename']; r['fasta_path']=str(p); r['sequence_check_status']='VALID' if p.exists() and fasta_valid(p) else 'NOT_VALIDATED';
  if r['reference_status']=='ACCEPTED' and r['sequence_check_status']=='VALID': valid.append(p)
  elif r['reference_status']=='ACCEPTED': issues.append({'reference_id':r['reference_id'],'issue':'missing_or_invalid_fasta','fasta_path':str(p)})
 with open(out/'reference_manifest.validated.tsv','w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]) if rows else ['reference_id'],delimiter='\t'); w.writeheader(); w.writerows(rows)
 if issues: json.dump(issues,open(out/'reference_validation_issues.json','w'),indent=2)
 if not valid: return {'status':'NOT_COMPLETED','reference_count':0,'issues':len(issues)}
 refs=out/'references.txt'; refs.write_text('\n'.join(str(p) for p in valid)+'\n')
 result={'status':'VALIDATED','reference_count':len(valid),'reference_list':str(refs.resolve())}
 if build_skani:
  exe=shutil.which('skani')
  if not exe: raise RuntimeError('skani not found on Linux server')
  version=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True).stdout.strip(); prefix=out/'skani_db'; subprocess.run(build_skani_command(exe,refs,prefix,threads),check=True); result.update({'status':'COMPLETED','skani_version':version,'index_path':str(prefix.resolve())})
 json.dump(result,open(out/'index_manifest.json','w'),indent=2); return result
