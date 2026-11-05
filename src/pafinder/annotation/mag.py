"""Phase 3A MAG protein annotation: Pyrodigal -> DIAMOND -> metadata joins."""
from pathlib import Path
import csv,json,hashlib,subprocess,shutil,re

def read_fasta(path):
 h=None; seq=[]
 with open(path) as f:
  for line in f:
   line=line.strip()
   if line.startswith('>'):
    if h is not None: yield h,''.join(seq)
    h=line[1:]; seq=[]
   else: seq.append(line)
  if h is not None: yield h,''.join(seq)

def predict_orfs(mag_fasta,outdir,mode='meta',threads=8):
 try: import pyrodigal
 except ImportError as e: raise RuntimeError('Pyrodigal is required; install it in the pafinder Conda environment') from e
 outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); faa=outdir/'genes.faa'; fna=outdir/'genes.fna'; gff=outdir/'genes.gff'; rows=[]
 with open(faa,'w') as ap, open(fna,'w') as an, open(gff,'w') as ag:
  ag.write('##gff-version 3\n')
  for contig,seq in read_fasta(mag_fasta):
   mag=Path(mag_fasta).stem; c=contig.split()[0]
   finder=pyrodigal.GeneFinder(meta=(mode=='meta'))
   genes = finder.find_genes(seq.encode())
   for i,g in enumerate(genes,1):
    gid=f'{mag}__{c}__gene_{i:05d}'; start,end,strand=g.begin,g.end,g.strand
    header=f'{gid} MAG={mag} contig={c} gene={gid}'; prot=str(g.translate()); nuc=seq[g.begin-1:g.end] if strand==1 else _revcomp(seq[g.begin-1:g.end])
    ap.write(f'>{header}\n{prot}\n'); an.write(f'>{header}\n{nuc}\n')
    attrs=f'ID={gid};MAG={mag};contig={c};gene_id={gid}'
    ag.write(f'{c}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{"+" if strand==1 else "-"}\t0\t{attrs}\n')
    rows.append({'mag_id':mag,'contig_id':c,'gene_id':gid,'start':start,'end':end,'strand':'+' if strand==1 else '-'})
 with open(outdir/'genes.tsv','w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=rows[0].keys() if rows else ['mag_id','contig_id','gene_id','start','end','strand'],delimiter='\t'); w.writeheader(); w.writerows(rows)
 return outdir/'genes.faa',outdir/'genes.fna',outdir/'genes.gff',outdir/'genes.tsv',len(rows)

def _revcomp(s): return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]
def sha256(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def build_diamond_command(exe, query, db, out, threads=8, **opts):
    """Build one DIAMOND blastp command with normalized options."""
    fields = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "bitscore"]
    cmd = [str(exe), "blastp", "--query", str(query), "--db", str(db), "--out", str(out), "--outfmt", "6", *fields, "--threads", str(threads)]
    sensitivity = opts.get("sensitivity")
    if sensitivity:
        flag = {"sensitive": "--sensitive", "very-sensitive": "--very-sensitive"}.get(str(sensitivity).lower())
        if not flag:
            raise ValueError(f"Unsupported DIAMOND sensitivity mode: {sensitivity}")
        cmd.append(flag)
    for key, value in opts.items():
        if value is None or key in {"sensitivity", "outfmt"}:
            continue
        cmd.extend([f"--{key.replace('_', '-')}", str(value)])
    return cmd

def diamond_search(query,ref,outdir,label,threads=8,**opts):
 exe=shutil.which('diamond')
 if not exe: raise RuntimeError('DIAMOND executable not found; install diamond from bioconda in the pafinder environment')
 outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); db=outdir/f'{label}.dmnd'; subprocess.run([exe,'makedb','--in',str(ref),'--db',str(db.with_suffix(''))],check=True)
 out=outdir/f'{label}_hits.raw.tsv'
 cmd=build_diamond_command(exe,query,db.with_suffix(''),out,threads,**opts)
 subprocess.run(cmd,check=True); return out

HIT_FIELDS = "qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore".split()
OUTPUT_HIT_FIELDS = HIT_FIELDS + ["query_span_coverage", "subject_span_coverage", "metadata_status"]

def join_hits(raw, metadata, outfile):
    """Join raw DIAMOND HSPs to metadata without dropping unmapped hits."""
    with open(metadata, newline="") as h:
        metadata_rows = list(csv.DictReader(h, delimiter="\t"))
    md = {r.get("record_id", ""): r for r in metadata_rows if r.get("record_id")}
    metadata_fields = sorted({k for r in metadata_rows for k in r if k != "record_id"})
    prefixed = [f"metadata_{k}" for k in metadata_fields]
    fields = OUTPUT_HIT_FIELDS + prefixed
    rows = []
    with open(raw) as h:
        for line in h:
            vals = line.rstrip("\n").split("\t")
            if len(vals) != len(HIT_FIELDS):
                continue
            d = dict(zip(HIT_FIELDS, vals))
            qlen, slen = max(1, int(d["qlen"])), max(1, int(d["slen"]))
            qspan = abs(int(d["qend"]) - int(d["qstart"])) + 1
            sspan = abs(int(d["send"]) - int(d["sstart"])) + 1
            d["query_span_coverage"] = f"{100*qspan/qlen:.3f}"
            d["subject_span_coverage"] = f"{100*sspan/slen:.3f}"
            match = md.get(d["sseqid"])
            d["metadata_status"] = "MAPPED" if match else "UNMAPPED"
            if match:
                d.update({f"metadata_{k}": v for k, v in match.items() if k != "record_id"})
            rows.append(d)
    with open(outfile, "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    return {"alignment_rows": len(rows), "output_hit_rows": len(rows), "metadata_mapped_rows": sum(r["metadata_status"] == "MAPPED" for r in rows), "unmapped_subject_ids": sorted({r["sseqid"] for r in rows if r["metadata_status"] == "UNMAPPED"})}

def write_empty_hits(path, metadata):
    """Create a valid, header-only hit table for zero ORFs or zero hits."""
    with open(metadata, newline="") as h:
        fields = sorted({k for r in csv.DictReader(h, delimiter="\t") for k in r if k != "record_id"})
    with open(path, "w", newline="") as o:
        csv.writer(o, delimiter="\t").writerow(OUTPUT_HIT_FIELDS + [f"metadata_{k}" for k in fields])

def run(config,mag):
 import yaml
 cfg=yaml.safe_load(open(config)); root=Path(__file__).resolve().parents[3]; out=(root/cfg.get('output_dir','build/phase3a')).resolve(); out.mkdir(parents=True,exist_ok=True)
 # Remove stale result artifacts before every run.
 for pattern in ('*_hits.raw.tsv', 'vfdb.dmnd', 'sarg.dmnd', 'vfdb', 'sarg'):
  for p in out.glob(pattern):
   if p.is_file(): p.unlink()
 faa,_,_,genes,predicted = predict_orfs(mag,out,cfg.get('mode','meta'),int(cfg.get('threads',8)))
 vf_meta=root/'build/databases/v1/vfdb/metadata.tsv'; arg_meta=root/'build/databases/v1/sarg_plus/metadata.tsv'
 if predicted == 0:
  for p,md in [(out/'vf_hits.tsv',vf_meta),(out/'arg_hits.tsv',arg_meta)]: write_empty_hits(p,md)
  manifest={'mag':str(Path(mag).resolve()),'mode':cfg.get('mode','meta'),'threads':cfg.get('threads',8),'predicted_genes':0,'vfdb_hit_count':0,'sarg_hit_count':0,'status':'NO_GENES_PREDICTED'}
  json.dump(manifest,open(out/'run_manifest.json','w'),indent=2); print('NO_GENES_PREDICTED: skipping DIAMOND'); return manifest
 raws={}; stats={}
 for label,ref in cfg['references'].items(): raws[label]=diamond_search(faa,root/ref,out,label,int(cfg.get('threads',8)),**cfg.get('diamond',{}))
 stats['vfdb']=join_hits(raws['vfdb'],vf_meta,out/'vf_hits.tsv'); stats['sarg']=join_hits(raws['sarg'],arg_meta,out/'arg_hits.tsv')
 manifest={'mag':str(Path(mag).resolve()),'mode':cfg.get('mode','meta'),'threads':cfg.get('threads',8),'predicted_genes':predicted,'vfdb_hit_count':stats['vfdb']['output_hit_rows'],'sarg_hit_count':stats['sarg']['output_hit_rows'],'metadata_mapping':stats,'status':'COMPLETED'}
 manifest['outputs']={p.name:sha256(p) for p in out.glob('*') if p.is_file()}; json.dump(manifest,open(out/'run_manifest.json','w'),indent=2); return manifest

