#!/usr/bin/env python3
import os,gzip,shutil,subprocess,time
GEN='./downloads/genomes'; SC='./scan/raw'; FAA='./scan/faa'
os.makedirs(SC,exist_ok=True); os.makedirs(FAA,exist_ok=True)
accs=sorted([l.strip() for l in open('./intermediate/accs45.txt') if l.strip()])
print("genomes:",len(accs),flush=True)
t0=time.time()
for i,acc in enumerate(accs):
    faa=os.path.join(FAA,acc+'.faa'); out=os.path.join(SC,acc+'.raw.tsv')
    if os.path.exists(out) and os.path.getsize(out)>0:
        print(f"[{i+1}] {acc} cached",flush=True); continue
    if not os.path.exists(faa):
        with gzip.open(os.path.join(GEN,acc,'protein.faa.gz'),'rb') as ig,open(faa,'wb') as og: shutil.copyfileobj(ig,og)
    r=subprocess.run(['python3','./scripts/scan_one.py',acc,faa,out,'24','1e-3'],capture_output=True,text=True)
    if r.returncode!=0:
        print(f"[{i+1}] {acc} FAIL {r.stderr[-300:]}",flush=True)
    else:
        print(f"[{i+1}] {acc} {r.stdout.strip()} ({time.time()-t0:.0f}s)",flush=True)
print("DONE",time.time()-t0,flush=True)
