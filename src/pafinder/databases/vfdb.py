import gzip,re,csv
from .common import fasta,write_fasta,stats
ID_RE=re.compile(r'(?P<vfg>VFG\d+).*?(?P<vf>VF\d+).*?(?P<vfc>VFC\d+)?',re.I)

def parse_ids(header):
 m=ID_RE.search(header); return ({'vfg_id':m.group('vfg').upper() if m else 'UNKNOWN','vf_id':m.group('vf').upper() if m else 'UNKNOWN','vfc_id':m.group('vfc').upper() if m and m.group('vfc') else 'UNKNOWN'})
def parse_xls(path):
 import xlrd
 with gzip.open(path,'rb') as h: book=xlrd.open_workbook(file_contents=h.read())
 sh=book.sheet_by_name('Bacterial VF (prototype)')
 hr=next(i for i in range(sh.nrows) if 'VFID' in [str(x).strip() for x in sh.row_values(i)])
 headers=[str(x).strip() for x in sh.row_values(hr)]
 rows=[]
 for i in range(hr+1,sh.nrows):
  vals=sh.row_values(i)
  if not any(str(x).strip() for x in vals): continue
  d=dict(zip(headers,vals)); d['xls_row_index_0based']=i; rows.append(d)
 return rows

def prepare(inp,out,evidence_dir=None):
 xls=parse_xls(inp/'VFDB/VFs.xls.gz'); byvf={str(r.get('VFID','')).strip().upper():r for r in xls}; rows=[]; records=[]; maps=[]
 for typ,rel,suf in [('protein','VFDB/VFDB_setA_pro.fas.gz','faa'),('nucleotide','VFDB/VFDB_setA_nt.fas.gz','fna')]:
  typerec=[]
  for i,r in enumerate(fasta(inp/rel)):
   ids=parse_ids(r.description); oid=r.id; rid=f'vfdb_{typ}_{i:06d}'; rec={'id':rid,'original_id':oid,'header':r.description,'sequence':str(r.seq),'sequence_type':typ,**ids}; typerec.append(rec); records.append(rec)
   md=byvf.get(ids['vf_id'],{}); rows.append({'db_name':'VFDB','record_id':rid,'original_id':oid,'original_header':r.description,'sequence_type':typ,**ids,'gene_name':'UNKNOWN','category':md.get('VFcategory','UNKNOWN'),'subcategory':'UNKNOWN','source_organism':md.get('Bacteria','UNKNOWN'),'evidence_provenance':'VFDB XLS','review_status':'UNREVIEWED','source_database':'VFDB Set A','database_version':'UNKNOWN','VFID':ids['vf_id'],'VF_name':md.get('VF_Name','UNKNOWN'),'Function':md.get('Function','UNKNOWN'),'Mechanism':md.get('Mechanism','UNKNOWN'),'Bacteria':md.get('Bacteria','UNKNOWN'),'VFCID':md.get('VFCID','UNKNOWN'),'Reference':md.get('Reference','UNKNOWN'),'xls_row_index_0based':md.get('xls_row_index_0based','UNKNOWN'),'xls_match_status':'MATCHED' if md else 'UNMATCHED'})
   write_fasta(typerec,out/f'vfdb/reference.{suf}')
 # evidence rules
 ev=[]
 if evidence_dir:
  p=evidence_dir/'vfdb_v1_evidence_rules.tsv'
  if p.exists():
   with open(p,encoding='utf-8') as h:
    for r in csv.DictReader(h,delimiter='\t'):
     if r.get('record_id')==r.get('record_id','').strip() and r.get('record_id')!='record_id': ev.append(r)
 evbyorig={r.get('original_id'):r for r in ev}; evbyrid={r.get('record_id'):r for r in ev}
 for r in rows:
  e=evbyorig.get(r['original_id']) or evbyrid.get(r['record_id'])
  if e:
   for k in ['evidence_tier_candidate','mechanism_role','module_dependency','evidence_rationale','review_priority','evidence_review_status','reference_url','review_status']:
    if e.get(k): r[k]=e[k]
  else: r.update({'evidence_tier_candidate':'UNREVIEWED','evidence_review_status':'UNREVIEWED','review_status':'UNREVIEWED'})
 statsall={'protein':stats([r for r in records if r['sequence_type']=='protein']),'nucleotide':stats([r for r in records if r['sequence_type']=='nucleotide']),'metadata_records':len(xls),'xls_vfid_count':len(byvf),'sequence_vfid_count':len(set(r['vf_id'] for r in records)),'metadata_mapped':sum(r['xls_match_status']=='MATCHED' for r in rows),'unmapped_sequences':sum(r['xls_match_status']!='MATCHED' for r in rows),'evidence_rules':len(ev),'evidence_matched':sum(1 for r in rows if r['sequence_type']=='protein' and r.get('evidence_tier_candidate') not in (None,'UNREVIEWED')),'vfc_conflicts':0}
 for r in rows: maps.append({k:r.get(k,'UNKNOWN') for k in ['record_id','original_id','sequence_type','vfg_id','vf_id','vfc_id','xls_row_index_0based','xls_match_status']})
 return rows,statsall,maps,ev
