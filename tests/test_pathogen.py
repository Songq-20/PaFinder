from pafinder.databases.pathogen import build_skani_command, search_skani_command

def test_skani_commands():
 assert build_skani_command('refs.txt','db',8)==['skani','sketch','-l','refs.txt','-o','db','-t','8']
 assert search_skani_command('db',['a.fna','b.fna'],'out',4)==['skani','search','-d','db','-t','4','-o','out','a.fna','b.fna']
