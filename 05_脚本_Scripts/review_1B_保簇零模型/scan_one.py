#!/usr/bin/env python3
"""Faithful replication of dbcan 5.2 HMM module (same as paper scan_proteome.py).
Z=875, included domains, i_evalue<domE(scan 1e-3), coverage=(hmm_to-hmm_from+1)/hmm_len>0.35."""
import sys, pyhmmer
HMM='./downloads/dbCAN.hmm'
def main():
    acc,faa,out,cpus=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4])
    domE=float(sys.argv[5]) if len(sys.argv)>5 else 1e-3
    covT=0.35
    hmms=[]
    with pyhmmer.plan7.HMMFile(HMM) as hf:
        for h in hf: hmms.append(h)
    Z=len(hmms); assert Z==875,f"profiles {Z}!=875"
    alphabet=pyhmmer.easel.Alphabet.amino()
    n=0
    with open(out,'w') as fo:
        fo.write("hmm_name\tfamily\thmm_length\ttarget_name\ttarget_length\ti_evalue\thmm_from\thmm_to\ttarget_from\ttarget_to\tcoverage\n")
        with pyhmmer.easel.SequenceFile(faa,digital=True,alphabet=alphabet) as seqs:
            targets=seqs.read_block()
            for hits in pyhmmer.hmmsearch(hmms,targets,cpus=cpus,domE=domE,Z=Z):
                for hit in hits:
                    for dom in hit.domains.included:
                        aln=dom.alignment
                        cov=(aln.hmm_to-aln.hmm_from+1)/aln.hmm_length
                        ie=dom.i_evalue
                        if ie<domE and cov>covT:
                            _hn=aln.hmm_name; name=_hn.decode() if isinstance(_hn,bytes) else _hn
                            fam=name[:-4] if name.endswith('.hmm') else name
                            fam='GT2' if fam.startswith('GT2_') else fam.split('_')[0]
                            _tn=aln.target_name; tn=_tn.decode() if isinstance(_tn,bytes) else _tn
                            fo.write(f"{name}\t{fam}\t{aln.hmm_length}\t{tn}\t{aln.target_length}\t{ie}\t{aln.hmm_from}\t{aln.hmm_to}\t{aln.target_from}\t{aln.target_to}\t{cov}\n")
                            n+=1
    print(f"{acc}\t{n}",flush=True)
main()
