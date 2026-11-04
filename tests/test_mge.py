from pafinder.annotation.mge import build_blastn_command

def test_blastn_command_uses_single_complete_outfmt():
    cmd = build_blastn_command('blastn', 'mag.fna', 'mgedb', 'hits.tsv', 8,
                               evalue='1e-5', max_target_seqs=50, task='megablast',
                               outfmt='bad duplicate')
    assert cmd.count('-outfmt') == 1
    i = cmd.index('-outfmt')
    assert cmd[i + 1] == '6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore'
    assert 'qseqid' not in cmd[i + 2:]
    assert '-max_target_seqs' in cmd
    assert '-max-target-seqs' not in cmd
    assert '-num_threads' in cmd and '-evalue' in cmd and '-task' in cmd
