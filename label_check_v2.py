#!/usr/bin/env python3
"""label_check_v2.py model.pdb [-v] [-ref REF.pdb]
-> 'labels=<pct> wrong=<n> maxD=<A> nbig=<n>' (plus the list of wrong pairs with -v).

R5 (claude_numbers, 27 Oct 2026). Differences from basin_study_2026-10-03/label_check.py:
  * reference = the TRUE ground truth (ground_truth.pdb next to this script = ground_factory5/onechain.pdb) (label_check.py used qfit_best/compare_to_best_holdfirst_slow/
    good.pdb, a phenix-2.0 refinement of the truth, up to 0.038 A off);
  * the A/B pair test treats X-ray-equivalent names as the same atom (each distance is the minimum over the atom's name and its
    partner's name), so a ring or carboxylate whose names are exchanged is not counted as a wrong pair;
  * His ND1/CD2 and CE1/NE2 and C-terminal O/OXT are added to the equivalences (JH 27 Oct: X-ray-equivalent solutions are not
    failures; the LGD driver searches the 4 flip combinations every cycle).
Heavy atoms only, waters excluded. A global A<->B relabel of the whole model is allowed: the smaller of (same, swapped) is reported
as wrong. maxD = for each model atom, the distance to the nearest truth atom of the same residue with the same or equivalent name
(any altloc); maxD is the largest such distance and nbig counts those over 1 A.
"""
import sys, math, os
# reference: ground_truth.pdb next to this script (UnTangle repo; identical to ground_factory5/onechain.pdb), else onechain.pdb
GT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ground_truth.pdb')
if not os.path.exists(GT): GT = '/home/jamesh/projects/build_challenge/1aho/ground_factory5/onechain.pdb'
args = sys.argv[1:]
verbose = '-v' in args
if '-ref' in args:
    GT = args[args.index('-ref') + 1]; del args[args.index('-ref'):args.index('-ref') + 2]
args = [a for a in args if a != '-v']
PAIRS = {'ASP': [('OD1', 'OD2')], 'GLU': [('OE1', 'OE2')], 'PHE': [('CD1', 'CD2'), ('CE1', 'CE2')], 'TYR': [('CD1', 'CD2'), ('CE1', 'CE2')],
         'ARG': [('NH1', 'NH2')], 'VAL': [('CG1', 'CG2')], 'LEU': [('CD1', 'CD2')], 'ASN': [('OD1', 'ND2')], 'GLN': [('OE1', 'NE2')],
         'HIS': [('ND1', 'CD2'), ('CE1', 'NE2')]}
EQ = {}
for r, ps in PAIRS.items():
    for a, b in ps: EQ[(r, a)] = b; EQ[(r, b)] = a
def eqname(resn, name, cterm):
    if name in ('O', 'OXT') and cterm: return 'OXT' if name == 'O' else 'O'
    return EQ.get((resn, name), '')
def read(fn):
    out = []
    for l in open(fn, errors='ignore'):
        if l.startswith('ENDMDL') and out: break
        if l[:6] not in ('ATOM  ', 'HETATM'): continue
        name = l[12:16].strip(); el = (l[76:78].strip() or name[0]); resn = l[17:20].strip()
        if el in ('H', 'D') or resn in ('HOH', 'WAT', 'DOD'): continue
        out.append((resn, int(l[22:26]), name, l[16], (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out
def d(a, b): return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))
gt = read(GT)
lastres = max(r for _, r, _, _, _ in gt)
gtab = {}; gti = {}
for resn, resi, name, alt, x in gt:
    gti.setdefault((resi, name), []).append(x)
    if alt in 'AB': gtab.setdefault((resi, name), {})[alt] = x
for fn in args:
    m = read(fn); mab = {}; worst = 0; worstat = ''; nbig = 0; alld = []
    for resn, resi, name, alt, x in m:
        if alt in 'AB': mab.setdefault((resi, name, resn), {})[alt] = x
        c = gti.get((resi, name), []) + gti.get((resi, eqname(resn, name, resi == lastres)), [])
        if c:
            dd = min(d(x, y) for y in c)
            if dd > worst: worst = dd; worstat = f'{resn}{resi}{alt.strip()}:{name}'
            alld.append((dd, f'{resn}{resi} {name} alt {alt.strip() or "_"}'))
            nbig += dd > 1.0
    same = []; swap = []
    for (resi, name, resn), v in mab.items():
        k = (resi, name)
        if not ('A' in v and 'B' in v and k in gtab and 'A' in gtab[k] and 'B' in gtab[k]): continue
        e = eqname(resn, name, resi == lastres)
        def cands(alt):
            c = [gtab[k][alt]]
            if e and (resi, e) in gtab and alt in gtab[(resi, e)]: c.append(gtab[(resi, e)][alt])
            return c
        dm = lambda x, alt: min(d(x, y) for y in cands(alt))
        ds = dm(v['A'], 'A') + dm(v['B'], 'B'); dw = dm(v['A'], 'B') + dm(v['B'], 'A')
        (same if ds <= dw else swap).append(f'{resn}{resi}:{name}')
    n = len(same) + len(swap)
    wrong = swap if len(swap) <= len(same) else same
    pct = max(len(same), len(swap)) / n * 100 if n else float('nan')
    pre = (fn + ' ') if len(args) > 1 else ''
    print(f"{pre}labels={pct:.1f} wrong={len(wrong)} maxD={worst:.2f} nbig={nbig} pairs={n} worst={worstat}")
    if verbose:
        alld.sort(reverse=True)
        print('  worst atoms: ' + '; '.join(f'{a} {v:.3f}' for v, a in alld[:3]))
    if verbose and wrong:
        res = {}
        for w in wrong:
            r, a = w.split(':'); res.setdefault(r, []).append(a)
        print('  wrong:', '; '.join(f"{r} {' '.join(a)}" for r, a in sorted(res.items(), key=lambda t: int(''.join(c for c in t[0] if c.isdigit())))))
