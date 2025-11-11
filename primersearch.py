###############################################################################
####                                                                       ####
####                      PRIMER CONCENSOUS FINDER                         ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################



from math import log2

def info_from_probs(p):
    H = 0
    for b in ('A','C','G','T'):
        if p[b] > 0:
            H -= p[b]*log2(p[b])
    return 2-H

def parse_hmm_probs(hmm_path):
    probs = []
    in_hmm = False
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith('HMM'):
                in_hmm = True
                next(fh)
                continue
            if not in_hmm:
                continue
            if line.strip() == '//' or not line.strip():
                break
            parts = line.split()
            if not parts[0].isdigit():
                continue
            pa, pc, pg, pt = map(float, parts[1:5])
            s = pa+pc+pg+pt
            if s == 0: p = {'A':.25,'C':.25,'G':.25,'T':.25}
            else: p = {b:v/s for b,v in zip('ACGT',(pa,pc,pg,pt))}
            probs.append(p)
    return probs

def best_window_within_consensus(hmm_path, fasta_path, w=80):

    probs = parse_hmm_probs(hmm_path)

    with open(fasta_path) as f:
        f.readline()
        consensus = f.read().replace('\n','').upper()
    L = min(len(probs), len(consensus))
    info = [info_from_probs(probs[i]) for i in range(L)]

    best = (-1,-1,-1)
    for i in range(L-w):
        mean_ic = sum(info[i:i+w])/w
        if mean_ic > best[2]:
            best = (i, i+w, mean_ic)
    start, end, mean_ic = best
    seq = consensus[start:end]
    return {'start': start+1, 'end': end, 'mean_IC_bits': round(mean_ic,3), 'sequence_nt': seq}

## Download Top10 DEGs from Dfam:


  
