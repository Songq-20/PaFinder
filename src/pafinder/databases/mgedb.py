import json
from .common import fasta,write_fasta,stats

def prepare(inp,out):
 d=inp/'MGEdb/mgedb/data'; db=json.load(open(d/'mge_records.json')); rows=[]; rec=[]; maps=[]
 for i,r in enumerate(fasta(d/'sequences.d/mge_records.fna')):
  rid=f'mgedb_nucleotide_{i:06d}'; rec.append({'id':rid,'original_id':r.id,'header':r.description,'sequence':str(r.seq),'sequence_type':'nucleotide'}); rows.append({'db_name':'MGEdb','record_id':rid,'original_id':r.id,'original_header':r.description,'sequence_type':'nucleotide','gene_name':'UNKNOWN','category':'UNKNOWN','subcategory':'UNKNOWN','source_organism':'UNKNOWN','evidence_provenance':'MGEdb','review_status':'UNKNOWN','source_database':'MGEdb','database_version':'UNKNOWN','Element_class':'UNKNOWN','Element_name':'UNKNOWN','Parent_element_ID':'UNKNOWN','Nomenclature':'UNKNOWN','Reference_information':'UNKNOWN'}); maps.append({'record_id':rid,'original_id':r.id,'sequence_type':'nucleotide'})
 write_fasta(rec,out/'mgedb/reference.fna')
 cds=[]
 for i,r in enumerate(fasta(d/'sequences.d/mge_cds.faa')):
  rid=f'mgedb_cds_{i:06d}'; cds.append({'id':rid,'original_id':r.id,'header':r.description,'sequence':str(r.seq),'sequence_type':'protein'}); rows.append({'db_name':'MGEdb','record_id':rid,'original_id':r.id,'original_header':r.description,'sequence_type':'protein','gene_name':'UNKNOWN','category':'CDS','subcategory':'UNKNOWN','source_organism':'UNKNOWN','evidence_provenance':'MGEdb','review_status':'UNKNOWN','source_database':'MGEdb','database_version':'UNKNOWN','Element_class':'UNKNOWN','Element_name':'UNKNOWN','Parent_element_ID':'UNKNOWN','Nomenclature':'UNKNOWN','Reference_information':'UNKNOWN'}); maps.append({'record_id':rid,'original_id':r.id,'sequence_type':'protein'})
 write_fasta(cds,out/'mgedb/cds_reference.faa')
 return rows,{'nucleotide':stats(rec),'cds':stats(cds),'json_records':len(db)},maps

def repair_metadata(inp, out):
 import json,csv
 from pathlib import Path
 d=inp/'MGEdb/mgedb/data'; records=json.load(open(d/'mge_records.json')); nom=json.load(open(d/'mge_nomenclature.json')); out=Path(out); rows=[]; maps=[]
 type_by_name={n:name for name,vals in nom.items() for n in vals}
 for i,r in enumerate(fasta(d/'sequences.d/mge_records.fna')):
  oid=r.id; element=oid.split('|')[0]; rec=records.get(element,{}); category=rec.get('type','NOT_PROVIDED'); name=rec.get('name') or element; nomen=type_by_name.get(name, name if name in records else 'NOT_PROVIDED'); accession=oid.split('|')[-1] if '|' in oid else 'NOT_PROVIDED'
  rid=f'mgedb_nucleotide_{i:06d}'; row={'db_name':'MGEdb','record_id':rid,'original_id':oid,'original_header':r.description,'sequence_type':'nucleotide','gene_name':'NOT_PROVIDED','category':category,'subcategory':rec.get('family','NOT_PROVIDED'),'source_organism':'; '.join(x.get('name','') for x in rec.get('organism',[])) or 'NOT_PROVIDED','evidence_provenance':rec.get('evidence_lvl','NOT_PROVIDED'),'review_status':'UNKNOWN','source_database':'MGEdb','database_version':'UNKNOWN','Element_class':category,'Element_name':name,'Parent_element_ID':rec.get('family','NOT_PROVIDED'),'Nomenclature':nomen,'Reference_information':';'.join(map(str,rec.get('references',[]))) or 'NOT_PROVIDED','reference_accession':accession,'metadata_status':'MAPPED' if rec else 'UNMAPPED'}; rows.append(row); maps.append({'record_id':rid,'original_id':oid,'sequence_type':'nucleotide','mge_record_id':element,'metadata_status':row['metadata_status']})
 def tab(rs,p):
  with open(p,'w',newline='') as h:w=csv.DictWriter(h,fieldnames=sorted({k for x in rs for k in x}),delimiter='\t'); w.writeheader(); w.writerows(rs)
 tab(rows,out/'metadata.tsv'); tab(maps,out/'id_mapping.tsv'); return rows
