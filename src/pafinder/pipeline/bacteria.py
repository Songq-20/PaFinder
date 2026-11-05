from pathlib import Path
import csv,json,shutil,subprocess
try:
 import yaml
except ImportError: yaml=None
from pafinder.annotation.mag import predict_orfs,diamond_search,join_hits
from pafinder.annotation.mge import blastn_search,parse_hits,proximity
from pafinder.integration.final import run as integrate
from pafinder.databases.pathogen import search_skani_command,parse_skani

def _files(inp):
 p=Path(inp); return [p] if p.is_file() else sorted(x for x in p.iterdir() if x.suffix.lower() in ('.fa','.fna','.fasta'))
def _db(db):
 d=Path(db); return {'vf':d/'build/databases/v1/vfdb/reference.faa','arg':d/'build/databases/v1/sarg_plus/reference.faa','mge':d/'build/databases/v1/mgedb/reference.fna','vfmeta':d/'build/databases/v1/vfdb/metadata.tsv','argmeta':d/'build/databases/v1/sarg_plus/metadata.tsv','mgemeta':d/'build/databases/v1/mgedb/metadata.tsv','pathogen':d/'build/databases/v1/pathogen_references'}
def run_pipeline(inp,out,db,threads=8,plot=False,skip=None,dry_run=False,skip_skani=False,config_root=None):
 out=Path(out); out.mkdir(parents=True,exist_ok=True); paths=_db(db); statuses=[]
 cfgroot=Path(config_root or Path(__file__).resolve().parents[3])/'config'
 p3a=yaml.safe_load(open(cfgroot/'phase3a.yaml')) if yaml else {'diamond':{'sensitivity':'sensitive','evalue':'1e-5','max_target_seqs':25}}
 p3b=yaml.safe_load(open(cfgroot/'phase3b.yaml')) if yaml else {'blastn':{'task':'megablast','evalue':'1e-5','max_target_seqs':50},'neighbor_window':5000}
 run_params={'diamond':p3a.get('diamond',{}),'blastn':p3b.get('blastn',{}),'neighbor_window':p3b.get('neighbor_window',5000)}
 seen_ids=set()
 for mag in _files(inp):
  mid=mag.stem
  if mid in seen_ids: raise ValueError(f'Input MAG ID collision: {mid}')
  seen_ids.add(mid)
  magroot=out/'mags'/mid; w=magroot/'work'; a=w/'phase3a'; b=w/'phase3b'; a.mkdir(parents=True,exist_ok=True); b.mkdir(parents=True,exist_ok=True)
  st={'mag_id':mid,'status':'FAILED','error':''}
  try:
   if dry_run: st.update(status='NOT_RUN'); statuses.append(st); continue
   faa,_,_,_,n=predict_orfs(mag,a,'meta',threads); st['predicted_genes']=n
   if n==0: st.update(status='SUCCESS_NO_HITS'); statuses.append(st); continue
   vfraw=diamond_search(faa,paths['vf'],a,'vfdb',threads,**run_params['diamond']); argr=diamond_search(faa,paths['arg'],a,'sarg',threads,**run_params['diamond'])
   join_hits(vfraw,paths['vfmeta'],a/'vf_hits.tsv'); join_hits(argr,paths['argmeta'],a/'arg_hits.tsv')
   mraw=blastn_search(mag,paths['mge'],b,threads,**run_params['blastn']); parse_hits(mraw,paths['mgemeta'],b/'mge_hits.tsv',{},mid)
   proximity(a/'genes.tsv',a/'vf_hits.tsv',b/'mge_hits.tsv',b/'vf_mge_proximity.tsv',run_params['neighbor_window']); proximity(a/'genes.tsv',a/'arg_hits.tsv',b/'mge_hits.tsv',b/'arg_mge_proximity.tsv',run_params['neighbor_window'])
   integrate(a,b,magroot,plot)
   if skip_skani: st['skani_status']='NOT_RUN'
   else:
    skdb=paths['pathogen']/'skani_db'; manifest=paths['pathogen']/'reference_manifest.tsv'
    if not skdb.exists() or not manifest.exists(): raise RuntimeError('skani index or reference_manifest.tsv missing; use --skip-skani to continue explicitly')
    raw=w/'skani_hits.raw.tsv'; exe=shutil.which('skani')
    if not exe: raise RuntimeError('skani executable missing; use --skip-skani to continue explicitly')
    subprocess.run(search_skani_command(skdb,[mag],raw,threads),check=True); matches=parse_skani(raw,manifest);
    ev=magroot/'evidence'; ev.mkdir(parents=True,exist_ok=True);
    with open(ev/'pathogen_matches.tsv','w',newline='') as fh:
     fields=list(matches[0]) if matches else ['query_mag_id','reference_id','ANI','AF_reference','AF_query','metadata_status']; ww=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); ww.writeheader(); ww.writerows(matches)
    st['skani_status']='SUCCESS_WITH_HITS' if matches else 'SUCCESS_NO_HITS'; st['pathogen_matches']=matches
   st.update(status='SUCCESS_WITH_HITS')
  except Exception as e: st['error']=str(e)
  statuses.append(st)
 Path(out/'qc').mkdir(exist_ok=True)
 # Aggregate one summary row per input MAG.
 allrows=[]
 for st in statuses:
  p=out/'mags'/st['mag_id']/'summary.tsv'
  if p.exists(): allrows.extend(csv.DictReader(open(p),delimiter='\t'))
 if allrows:
  fields=list(allrows[0]); sk={x['mag_id']:x.get('pathogen_matches',[]) for x in statuses}
  for row in allrows:
   ms=sk.get(row['mag_id'],[]); row['pathogen_match_count']=len(ms); row['closest_reference']=ms[0].get('reference_id','') if ms else ''; row['closest_ANI']=ms[0].get('ANI','') if ms else ''; row['closest_AF_reference']=ms[0].get('AF_reference','') if ms else ''; row['closest_AF_query']=ms[0].get('AF_query','') if ms else ''
   for k in ['pathogen_match_count','closest_reference','closest_ANI','closest_AF_reference','closest_AF_query']:
    if k not in fields: fields.append(k)
  with open(out/'summary.tsv','w',newline='') as fh: ww=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); ww.writeheader(); ww.writerows(allrows)
  with open(out/'summary_stats.tsv','w',newline='') as fh: ww=csv.writer(fh,delimiter='\t'); ww.writerow(['metric','value']); ww.writerow(['mag_count',len(allrows)])
 with open(out/'qc'/'mag_status.tsv','w',newline='') as h: wri=csv.DictWriter(h,fieldnames=['mag_id','status','predicted_genes','error'],delimiter='\t'); wri.writeheader(); wri.writerows({k:r.get(k,'') for k in wri.fieldnames} for r in statuses)
 json.dump({'input':str(Path(inp).resolve()),'database':str(Path(db).resolve()),'threads':threads,'parameters':run_params,'mag_count':len(statuses),'statuses':statuses},open(out/'qc'/'run_manifest.json','w'),indent=2)
 return statuses
