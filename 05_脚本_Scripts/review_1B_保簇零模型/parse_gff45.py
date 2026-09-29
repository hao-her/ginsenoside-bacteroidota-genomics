#!/usr/bin/env python3
"""Parse RefSeq genomic.gff.gz -> gene order table (same as paper parse_gff.py)."""
import gzip, os
GEN='./downloads/genomes'; OUTD='./scan/gff'
os.makedirs(OUTD,exist_ok=True)
def parse(acc):
    out=os.path.join(OUTD,acc+'.genes.tsv')
    rows=[]; seqorder={}
    with gzip.open(os.path.join(GEN,acc,'genomic.gff.gz'),'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9 or p[2]!='CDS': continue
            attrs=dict(kv.split('=',1) for kv in p[8].split(';') if '=' in kv)
            pid=attrs.get('protein_id','')
            if not pid: continue
            if p[0] not in seqorder: seqorder[p[0]]=len(seqorder)
            rows.append((seqorder[p[0]],p[0],int(p[3]),int(p[4]),p[6],attrs.get('locus_tag',''),pid))
    rows.sort(key=lambda r:(r[0],r[2]))
    with open(out,'w') as fo:
        fo.write("gene_ordinal\treplicon\tstart\tend\tstrand\tlocus_tag\tprotein_id\n")
        for i,r in enumerate(rows):
            fo.write(f"{i}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\n")
    return len(rows)
accs=sorted([l.strip() for l in open('./intermediate/accs45.txt') if l.strip()])
tot=0
for acc in accs:
    n=parse(acc); tot+=n
print(f"parsed {len(accs)} genomes, {tot} CDS total")
