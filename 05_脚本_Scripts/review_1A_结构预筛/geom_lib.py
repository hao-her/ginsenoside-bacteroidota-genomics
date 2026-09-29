"""Catalytic-pair geometry library for structure pre-screen (task 1A).
Maps the two structure-anchored catalytic Glu (GH3 cols 2165/3114 grounded on
BlBG3 5Z9S Glu162/Glu524; GH5 cols 1241/1473 grounded on BcelFp 6KDD Glu144/Glu260)
onto each candidate by pairwise alignment to the structural reference, then measures
the acid/base--nucleophile carboxyl geometry in the predicted model."""
import json, numpy as np, warnings
warnings.simplefilter('ignore')
from Bio.PDB import PDBParser, MMCIFParser
from Bio.SeqUtils import seq1
from Bio import Align
from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")
def make_aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = BLOSUM62
    a.open_gap_score = -11; a.extend_gap_score = -1
    a.mode = 'global'
    return a
ALN = make_aligner()

def transfer_positions(ref_seq, ref_idx_list, cand_seq):
    """Return candidate 0-based indices aligned to given ref 0-based indices."""
    aln = ALN.align(ref_seq, cand_seq)[0]
    # aln.aligned: blocks [[(rs,re),...],[(cs,ce),...]]
    ref_blocks, cand_blocks = aln.aligned
    # build ref_index -> cand_index map
    m = {}
    for (rs,re),(cs,ce) in zip(ref_blocks, cand_blocks):
        for k in range(re-rs):
            m[rs+k] = cs+k
    out = [m.get(i, None) for i in ref_idx_list]
    score = aln.score
    return out, score

def parse_struct(path):
    if path.endswith('.cif'):
        p = MMCIFParser(QUIET=True)
    else:
        p = PDBParser(QUIET=True)
    return p.get_structure('s', path)[0]

def chain_residues(model, chain_id=None):
    """Return list of standard residues in order, plus seq string."""
    ch = model[chain_id] if chain_id else list(model)[0]
    res = [r for r in ch if r.id[0]==' ']
    seq = ''.join(seq1(r.get_resname()) for r in res)
    return res, seq

def carboxyl_O(res):
    names = {'GLU':['OE1','OE2'],'ASP':['OD1','OD2']}
    return [res[a] for a in names.get(res.get_resname(),[]) if a in res]
def cd_atom(res):
    n = {'GLU':'CD','ASP':'CG'}.get(res.get_resname())
    return res[n] if (n and n in res) else None
def min_OO(r1,r2):
    a1,a2=carboxyl_O(r1),carboxyl_O(r2)
    return float(min(np.linalg.norm(x.coord-y.coord) for x in a1 for y in a2)) if a1 and a2 else None
def cd_cd(r1,r2):
    c1,c2=cd_atom(r1),cd_atom(r2)
    return float(np.linalg.norm(c1.coord-c2.coord)) if c1 is not None and c2 is not None else None
def ca_ca(r1,r2):
    if 'CA' in r1 and 'CA' in r2:
        return float(np.linalg.norm(r1['CA'].coord-r2['CA'].coord))
    return None
def plddt(res):
    # AF stores per-atom pLDDT in B-factor; use CA
    if 'CA' in res: return float(res['CA'].bfactor)
    return float(np.mean([a.bfactor for a in res]))

# reference gold-standard distances (measured from PDB)
REF_GEOM = {
 'GH3': {'pdb':'5Z9S','ab':162,'nu':524,'minOO':6.26,'cdcd':7.39,'caca':10.96},
 'GH5': {'pdb':'6KDD','ab':144,'nu':260,'minOO':3.35,'cdcd':4.80,'caca':10.17},
}
