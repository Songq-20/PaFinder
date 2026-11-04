import gzip,hashlib,json,re
from pathlib import Path

def fasta(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',encoding='utf-8',errors='replace') as h:
        header=None; seq=[]
        for line in h:
            line=line.rstrip('\n')
            if line.startswith('>'):
                if header is not None: yield type('Rec',(),{'id':header.split()[0],'description':header,'seq':''.join(seq)})()
                header=line[1:]; seq=[]
            else: seq.append(line.strip())
        if header is not None: yield type('Rec',(),{'id':header.split()[0],'description':header,'seq':''.join(seq)})()
def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def write_fasta(records,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w') as h:
        for r in records: h.write(f'>{r["id"]} {r.get("header","")}\n{r["sequence"]}\n')
def stats(records):
    seq=[r['sequence'] for r in records]; ids=[r['original_id'] for r in records]
    return {'records':len(seq),'valid_sequences':sum(bool(re.fullmatch(r'[A-Za-z*.-]+',s)) for s in seq),'unique_ids':len(set(ids)),'duplicate_ids':len(ids)-len(set(ids)),'identical_sequences':len(seq)-len(set(seq)),'length_min':min(map(len,seq),default=0),'length_max':max(map(len,seq),default=0)}
