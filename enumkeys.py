#!/usr/bin/env python3
"""
Enumerate orbit-key equivalence classes for brobdicube reduction steps,
then generate costalgs files for each useful face-subset key.

Orbit key: 6 characters, one per face in U L F R B D order.
  digit  -- color group (0, 1, 2, ..., up to 5)
  'x'    -- already solved; algorithms must not disturb this face

Face-subset keys  {0,x}^6     : reduced by the 24 cube rotations only.
Reduction keys    {0..5,x}^6  : also reduced by all digit relabelings
                                 (i.e., permutations of the digit labels).

State space size for a reduction key = n_act! / prod(n_g!)
  where n_act = 4 * (active faces) and n_g = 4 * (faces in group g),
  since every center orbit has 4 sticker positions per face and all
  permutations of active stickers are reachable by center 3-cycles.
"""

import os
import sys
from math import factorial
from collections import Counter

# Face indices: U=0 L=1 F=2 R=3 B=4 D=5

# Rotation permutations: new_key[i] = old_key[perm[i]]
# Three 90-degree generators (clockwise as seen from the indicated face):
X = (4, 1, 0, 3, 5, 2)  # R face: U->F->D->B->U  (L, R fixed)
Y = (0, 4, 1, 2, 3, 5)  # U face: F->R->B->L->F  (U, D fixed)
Z = (1, 5, 2, 0, 4, 3)  # F face: U->R->D->L->U  (F, B fixed)

def compose(a, b):
    return tuple(a[b[i]] for i in range(6))

# Generate all 24 rotation group elements via BFS
rots = [tuple(range(6))]
seen = {tuple(range(6))}
queue = [tuple(range(6))]
while queue:
    p = queue.pop(0)
    for g in (X, Y, Z):
        np = compose(g, p)
        if np not in seen:
            seen.add(np)
            rots.append(np)
            queue.append(np)

assert len(rots) == 24

def apply_rot(perm, key):
    return ''.join(key[perm[i]] for i in range(6))

def relabel(key):
    """Canonicalize digit labels by order of first appearance (lex-min under
    all digit permutations): the first digit seen becomes '0', the second
    distinct digit becomes '1', etc.  'x' is left unchanged."""
    mapping = {}
    next_d = 0
    result = []
    for c in key:
        if c == 'x':
            result.append('x')
        else:
            if c not in mapping:
                mapping[c] = str(next_d)
                next_d += 1
            result.append(mapping[c])
    return ''.join(result)

def canonical_subset(key):
    """Canonical form for face-subset keys: lex-min over rotations."""
    return min(apply_rot(r, key) for r in rots)

def canonical_reduction(key):
    """Canonical form for reduction keys: lex-min over rotations whose
    raw face-subset exactly equals the canonical face-subset, with digit
    relabeling.  This guarantees every canonical reduction key's 'x'
    positions match those of the canonical face-subset key for that orbit,
    so the input costalgs file is always expanded-{canon_fs}.costalgs."""
    raw_fs   = ''.join('0' if c != 'x' else 'x' for c in key)
    canon_fs = canonical_subset(raw_fs)
    cands = [relabel(apply_rot(r, key))
             for r in rots
             if ''.join('0' if c != 'x' else 'x' for c in apply_rot(r, key)) == canon_fs]
    return min(cands)

def space_size(key):
    """Multinomial coefficient: n_act! / prod(n_g!) with 4 positions per face."""
    counts = Counter(c for c in key if c != 'x')
    n_act = sum(counts.values()) * 4
    result = factorial(n_act)
    for face_count in counts.values():
        result //= factorial(face_count * 4)
    return result

# ── 1. Face-subset keys {0, x}^6 ────────────────────────────────────────────

all_face_keys = set()
for mask in range(64):   # 2^6 = 64
    key = ''.join('0' if (mask >> i) & 1 else 'x' for i in range(6))
    all_face_keys.add(canonical_subset(key))

# Exclude trivial: xxxxxx (no active faces) and 0xxxxx (1 active face, nothing to do)
SKIP_FACE = {'xxxxxx', '0xxxxx'}
face_keys = sorted(k for k in all_face_keys if k not in SKIP_FACE)

print("Face-subset keys (for algutil -f; '0'=active, 'x'=fixed):")
for k in face_keys:
    n = k.count('0')
    print(f"  {k}  ({n} active face{'s' if n != 1 else ''})")

print()

# ── 2. Reduction keys {0..5, x}^6 ───────────────────────────────────────────
# 7 values per position (x, 0, 1, 2, 3, 4, 5), 7^6 = 117649 total.
# After relabeling, digits are always 0..k-1 with k <= 6.

CHARS = 'x012345'
red_keys = set()
for n in range(7**6):
    tmp, key = n, ''
    for _ in range(6):
        key += CHARS[tmp % 7]
        tmp //= 7
    red_keys.add(canonical_reduction(key))

print("Reduction keys ('0','1',...=color groups, 'x'=fixed), sorted by state space size:")
for k in sorted(red_keys, key=lambda k: (space_size(k), k)):
    ss = space_size(k)
    if ss <= 1 or ss > 1_000_000_000:
        continue
    nact = sum(1 for c in k if c != 'x')
    ng   = len(set(c for c in k if c != 'x'))
    print(f"  {k}  space={ss:>20,}  ({nact} active, {ng} group{'s' if ng != 1 else ''})")

print()

# ── 3. Pipeline: generate costalgs files for each face-subset key ────────────

root = os.path.dirname(os.path.abspath(__file__))
os.chdir(root)
if os.system("make") != 0:
    sys.exit("make failed")

pipe_dir = os.path.join(root, "pipeline")
os.makedirs(pipe_dir, exist_ok=True)
os.chdir(pipe_dir)

A = "../algutil"

import time

print("Generating expanded.algs ...")
os.system(f"{A} -g -U < ../seed.algs 2>/dev/null | {A} -xu > expanded.algs 2>/dev/null")
with open("expanded.algs") as f:
    total = sum(1 for _ in f)
print(f"  {total:,} algorithms")
print()

F = "../findpermnet"

def face_subset_of(red_key):
    return ''.join('0' if c != 'x' else 'x' for c in red_key)

filtered = sorted(
    (k for k in red_keys if 1 < space_size(k) <= 1_000_000_000),
    key=lambda k: (space_size(k), k)
)

print("Generating costalgs files for each face-subset key:")
costalgs_counts = {}
for fsk in face_keys:
    if fsk == '000000':
        os.system(f"{A} -c -d 1000 -pc < expanded.algs > expanded-{fsk}.costalgs 2>/dev/null")
    else:
        os.system(f"{A} -f {fsk} < expanded.algs > expanded-{fsk}.algs 2>/dev/null")
        os.system(f"{A} -c -d 1000 -pc < expanded-{fsk}.algs > expanded-{fsk}.costalgs 2>/dev/null")
    with open(f"expanded-{fsk}.costalgs") as f:
        n = sum(1 for _ in f)
    costalgs_counts[fsk] = n
    print(f"  expanded-{fsk}.costalgs: {n:,} algorithms")

print()

def parse_permnet_log(logfile):
    """Return dict with: found, depth, cost, best_covered, total."""
    st = dict(found=False, depth=None, cost=None, best_covered=0, total=None)
    try:
        lines = open(logfile).readlines()
    except FileNotFoundError:
        return st
    in_found = False
    for line in lines:
        s = line.strip()
        if s.startswith('depth ') and ': best ' in s:
            try:
                frac = s.split(': best ')[1].split()[0]
                m, t = frac.split('/')
                st['best_covered'] = int(m)
                st['total'] = int(t)
            except (ValueError, IndexError):
                pass
        elif s.startswith('FOUND depth '):
            st['found'] = True
            in_found = True
            try:
                st['depth'] = int(s.split()[-1])
            except ValueError:
                pass
        elif in_found and s.startswith('cost ') and st['cost'] is None:
            try:
                st['cost'] = int(s.split()[1])
            except ValueError:
                pass
            in_found = False
    return st

max_ss = max(space_size(k) for k in filtered)

print("Running findpermnet for each reduction key (sorted by space size):")
for k in filtered:
    ss      = space_size(k)
    nact    = sum(1 for c in k if c != 'x')
    ng      = len(set(c for c in k if c != 'x'))
    fsk     = face_subset_of(k)
    ncands  = costalgs_counts.get(fsk, 0)
    infile  = f"expanded-{fsk}.costalgs"
    logfile = f"permnet-{k}.log"
    # Target: b * w * ss = 5 * 5 * max_ss (equal CPU budget per key)
    budget  = 25 * max_ss
    b = max(5, round(5 * (max_ss / ss) ** 0.5))
    w = b
    if ncands and b > ncands:
        b = ncands
        w = max(b, round(budget / (b * ss)))
    print(f"\n{k}  space={ss:,}  ({nact} active, {ng} group{'s' if ng!=1 else ''})"
          f"  b={b} w={w}  candidates={ncands:,}")
    sys.stdout.flush()
    t0 = time.time()
    os.system(f"{F} -b {b} -w {w} -i {k} < {infile} > {logfile} 2>&1")
    elapsed = time.time() - t0
    st = parse_permnet_log(logfile)
    if st['found']:
        print(f"  FOUND  depth={st['depth']}  cost={st['cost']:,}  ({elapsed:.1f}s)")
    else:
        cov = st['best_covered']
        tot = st['total'] or ss
        print(f"  NOT FOUND  best={cov:,}/{tot:,}  ({elapsed:.1f}s)")
    sys.stdout.flush()

# ── 4. Run pareto on all log files ──────────────────────────────────────────

P = "../pareto"
algs_dir = "algs3"

log_files = [f"permnet-{k}.log" for k in filtered]
found_logs = [lf for lf in log_files if parse_permnet_log(lf)['found']]
print(f"\nRunning pareto on {len(log_files)} log files"
      f" ({len(found_logs)} with solutions) -> {algs_dir}/")
sys.stdout.flush()
os.system(f"{P} -o {algs_dir} {' '.join(log_files)}")
n_algs = sum(1 for f in os.listdir(algs_dir) if f.endswith('.alg'))
print(f"  {algs_dir}/: {n_algs} .alg file{'s' if n_algs != 1 else ''}")

# ── 5. Brobdicube -e to precompute probability caches ───────────────────────

B = "../brobdicube"
CUBE_SIZE = 16384

print(f"\nPrecomputing probability caches: brobdicube -e -D {algs_dir} {CUBE_SIZE}")
sys.stdout.flush()
t0 = time.time()
os.system(f"{B} -e -D {algs_dir} {CUBE_SIZE}")
elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")

os.chdir(root)
