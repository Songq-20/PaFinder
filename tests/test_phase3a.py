from pathlib import Path
import pytest
from pafinder.annotation.mag import read_fasta

def test_small_fasta_reader():
 p=Path(__file__).parent/'fixtures'/'small_mag.fna'; assert list(read_fasta(p))==[('contig1','ATGAAATAGATGCCCCCC')]

def test_missing_pyrodigal_is_explicit(tmp_path,monkeypatch):
 import sys
 monkeypatch.setitem(sys.modules,'pyrodigal',None)
 from pafinder.annotation.mag import predict_orfs
 with pytest.raises(RuntimeError,match='Pyrodigal'):
  predict_orfs(Path(__file__).parent/'fixtures'/'small_mag.fna',tmp_path)

def test_diamond_command_normalizes_options():
 from pafinder.annotation.mag import build_diamond_command
 cmd=build_diamond_command('diamond','genes.faa','vfdb','hits.tsv',8,sensitivity='sensitive',evalue='1e-5',max_target_seqs=25,outfmt='bad')
 assert '--sensitive' in cmd and 'sensitivity' not in cmd
 assert cmd.count('--outfmt') == 1
 i=cmd.index('--outfmt'); assert cmd[i+1:i+4] == ['6','qseqid','sseqid']
 assert '--outfmt' not in cmd[i+1+1+12:]
 assert '--evalue' in cmd and '--max-target-seqs' in cmd

def test_diamond_very_sensitive_flag():
 from pafinder.annotation.mag import build_diamond_command
 cmd=build_diamond_command('diamond','q','db','o',sensitivity='very-sensitive')
 assert '--very-sensitive' in cmd

def test_join_hits_prefixes_metadata_and_uses_span_coverage(tmp_path):
 from pafinder.annotation.mag import join_hits
 md=tmp_path/'md.tsv'; raw=tmp_path/'raw.tsv'; out=tmp_path/'out.tsv'
 md.write_text('record_id\tgene_name\tevalue\nsub1\tfoo\tmeta\n')
 raw.write_text('q1\tsub1\t99\t8\t10\t17\t20\t2\t9\t40\t1e-5\t50\nq1\tmissing\t90\t5\t1\t5\t10\t1\t5\t10\t1\t20\n')
 stats=join_hits(raw,md,out); rows=list(__import__('csv').DictReader(open(out),delimiter='\t'))
 assert stats['alignment_rows']==2 and rows[0]['metadata_status']=='MAPPED' and rows[1]['metadata_status']=='UNMAPPED'
 assert rows[0]['metadata_evalue']=='meta' and rows[0]['query_span_coverage']=='40.000'

def test_empty_hit_table_has_header(tmp_path):
 from pafinder.annotation.mag import write_empty_hits
 md=tmp_path/'md.tsv'; out=tmp_path/'hits.tsv'; md.write_text('record_id\tgene_name\n')
 write_empty_hits(out,md); assert len(out.read_text().splitlines())==1

def test_mge_union_and_proximity(tmp_path):
 from pafinder.annotation.mge import parse_hits, proximity
 raw=tmp_path/'r'; md=tmp_path/'m'; out=tmp_path/'o'; genes=tmp_path/'g'; prox=tmp_path/'p'
 raw.write_text('contig1\tmgedb_nucleotide_000001\t99\t10\t1\t10\t100\t1\t10\t100\t1e-5\t20\ncontig1\tmgedb_nucleotide_000001\t99\t10\t11\t20\t100\t8\t17\t100\t1e-5\t20\n')
 md.write_text('record_id\tElement_name\n mgedb_nucleotide_000001\tX\n'.replace(' ','',1))
 rows,s=parse_hits(raw,md,out,{'partial_coverage':.1},'MAG')
 assert s[0]['merged_subject_covered_bp']=='20' and s[0]['hsp_count']=='2'
 genes.write_text('mag_id\tcontig_id\tgene_id\tstart\tend\nMAG\tcontig1\tg1\t15\t25\n'); fh=tmp_path/'fh'; fh.write_text('qseqid\n'+'g1\n'); assert proximity(genes,fh,out,prox,5)==1
