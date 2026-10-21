#!/usr/bin/env python3
"""
Score AF2-multimer GPCR-Galpha complexes: ipTM/pTM from score JSONs plus
interface metrics (contacts, interface pLDDT/PAE, pDockQ, alpha5 engagement)
computed from the top-ranked PDB.

Usage: python3 score_complexes.py <colabfold_output_dir> <out.tsv> [galpha.fasta] [gpcr.fasta]
"""
import json, glob, os, re, sys
import numpy as np
from math import log, exp

DIST_CUT, PLDDT_CUT, PAE_CUT = 8.0, 50.0, 15.0
ALPHA5_LEN = 26   # C-terminal helix window on Galpha

def read_fasta(p):
    s, n = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'): n = line[1:].split()[0]; s[n] = ''
        elif n: s[n] += line.strip()
    return s

def parse_pdb(pdb):
    ch_atoms = {}
    for line in open(pdb):
        if not line.startswith('ATOM'): continue
        ch = line[21]; resnum = int(line[22:26]); name = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        b = float(line[60:66])
        d = ch_atoms.setdefault(ch, {}).setdefault(resnum, {'coords': [], 'b': None})
        d['coords'].append(xyz)
        if name == 'CA' or d['b'] is None: d['b'] = b
    return ch_atoms

def pdockq_from(n_contacts, avg_plddt):
    if n_contacts == 0: return 0.0
    x = avg_plddt * log(n_contacts)
    return 0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018

def interface_metrics(pdb, pae_json):
    ch = parse_pdb(pdb)
    chains = sorted(ch)
    if len(chains) < 2: return None
    A, B = chains[0], chains[1]
    pae = np.array(json.load(open(pae_json))['predicted_aligned_error'], float)
    ra, rb = ch[A], ch[B]
    nA = max(ra)
    ca = {r: np.array(v['coords']) for r, v in ra.items()}
    cb = {r: np.array(v['coords']) for r, v in rb.items()}
    keysA, keysB = list(ca), list(cb)
    MA = np.array([ca[r].mean(0) for r in keysA])
    MB = np.array([cb[r].mean(0) for r in keysB])
    D = np.linalg.norm(MA[:, None, :] - MB[None, :, :], axis=-1)
    contacts = []
    for i, j in np.argwhere(D < DIST_CUT + 12):
        r1, r2 = keysA[i], keysB[j]
        d = np.linalg.norm(ca[r1][:, None, :] - cb[r2][None, :, :], axis=-1).min()
        if d > DIST_CUT: continue
        b1, b2 = ra[r1]['b'], rb[r2]['b']
        if b1 < PLDDT_CUT or b2 < PLDDT_CUT: continue
        gi, gj = r1 - 1, nA + r2 - 1
        p = min(pae[gi, gj], pae[gj, gi])
        if p > PAE_CUT: continue
        contacts.append((r1, r2, d, b1, b2, p))
    Lb = max(rb)
    if not contacts:
        return dict(n_contacts=0, n_res_gpcr=0, n_res_galpha=0,
                    avg_iface_plddt=np.nan, avg_iface_pae=np.nan, pdockq=0.0,
                    n_alpha5_contacts=0, frac_iface_alpha5=np.nan,
                    alpha5_res_engaged=0, min_iface_pae=np.nan,
                    iface_res_galpha='', iface_res_gpcr='')
    resA = sorted({c[0] for c in contacts}); resB = sorted({c[1] for c in contacts})
    plddts = [c[3] for c in contacts] + [c[4] for c in contacts]
    avg_plddt = float(np.mean(plddts))
    paes = [c[5] for c in contacts]
    a5_start = Lb - ALPHA5_LEN + 1
    a5_contacts = [c for c in contacts if c[1] >= a5_start]
    a5_res = sorted({c[1] for c in a5_contacts})
    return dict(
        n_contacts=len(contacts), n_res_gpcr=len(resA), n_res_galpha=len(resB),
        avg_iface_plddt=avg_plddt, avg_iface_pae=float(np.mean(paes)),
        min_iface_pae=float(np.min(paes)),
        pdockq=pdockq_from(len(contacts), avg_plddt),
        n_alpha5_contacts=len(a5_contacts),
        frac_iface_alpha5=len(a5_contacts) / len(contacts),
        alpha5_res_engaged=len(a5_res),
        iface_res_galpha=';'.join(map(str, resB)),
        iface_res_gpcr=';'.join(map(str, resA)))

def main(indir, outtsv):
    jobs = [os.path.basename(f).replace('.done.txt', '')
            for f in sorted(glob.glob(os.path.join(indir, '*.done.txt')))]
    rows = []
    for k, job in enumerate(jobs, 1):
        scj = sorted(glob.glob(os.path.join(indir, f'{job}_scores_rank_001*.json')))
        pdbs = sorted(glob.glob(os.path.join(indir, f'{job}_unrelaxed_rank_001*.pdb')))
        paej = os.path.join(indir, f'{job}_predicted_aligned_error_v1.json')
        if not (scj and pdbs and os.path.exists(paej)):
            print(f'[skip] {job}: missing outputs', flush=True); continue
        sc = json.load(open(scj[0]))
        model = re.search(r'model_(\d+)', os.path.basename(scj[0])).group(1)
        row = dict(job=job, top_model=int(model), iptm=sc['iptm'], ptm=sc['ptm'],
                   mean_plddt=float(np.mean(sc['plddt'])), max_pae=sc['max_pae'])
        # all-model ipTM spread
        allsc = sorted(glob.glob(os.path.join(indir, f'{job}_scores_rank_*.json')))
        iptms = [json.load(open(f))['iptm'] for f in allsc]
        row['iptm_max'] = max(iptms); row['iptm_min'] = min(iptms)
        row['iptm_mean'] = float(np.mean(iptms)); row['n_models'] = len(iptms)
        m = interface_metrics(pdbs[0], paej)
        if m: row.update(m)
        rows.append(row)
        if k % 20 == 0: print(f'  scored {k}/{len(jobs)}', flush=True)
    import csv
    if not rows: sys.exit('no scored complexes')
    cols = list(rows[0])
    with open(outtsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t'); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'wrote {len(rows)} rows -> {outtsv}')

if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
