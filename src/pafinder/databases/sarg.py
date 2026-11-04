from .common import fasta,write_fasta,stats
import csv

def prepare(inp,out):
 d=inp/'SARG+'; allr=[]; rows=[]; maprows=[]
 for source,fn in [('reference','sarg_ref.fa'),('extension','sarg_ext.fa')]:
  for i,r in enumerate(fasta(d/fn)):
   oid=r.id; rid=f'sarg_{source}_{i:06d}'; allr.append({'id':rid,'original_id':oid,'header':r.description,'sequence':str(r.seq),'sequence_type':'protein','source':source})
   parts=r.description.split('|'); rows.append({'db_name':'SARG+','record_id':rid,'original_id':oid,'original_header':r.description,'sequence_type':'protein','gene_name':parts[2] if len(parts)>2 else 'UNKNOWN','category':parts[1] if len(parts)>1 else 'UNKNOWN','subcategory':'UNKNOWN','source_organism':'UNKNOWN','evidence_provenance':source,'review_status':'UNKNOWN','source_database':'SARG+','database_version':'UNKNOWN','ARG_type':parts[1] if len(parts)>1 else 'UNKNOWN','ARG_subtype':parts[2] if len(parts)>2 else 'UNKNOWN','provenance':source})
   maprows.append({'record_id':rid,'original_id':oid,'source':source})
 write_fasta(allr,out/'sarg_plus/reference.faa')
 return rows,{'all':stats(allr),'reference_count':sum(x['source']=='reference' for x in allr),'extension_count':sum(x['source']=='extension' for x in allr),'metadata_records':0},maprows
