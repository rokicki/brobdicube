#!/usr/bin/env python3
"""
Enumerate orbit-key equivalence classes for brobdicube reduction steps,
then run the full reduction pipeline: costalgs generation, findpermnet beam
search, pareto extraction, and brobdicube probability-cache precomputation.

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

Usage:
    python3 runallreductions.py [options]

Options:
    --dir DIR        working directory for all pipeline output files
                     (default: pipeline, relative to this script)
    --max-cost N     filter bigseed.algs to cost <= N before expansion;
                     conjugations of the surviving seeds may still exceed N
                     and are kept (the filter is not applied to costalgs output)
"""

import argparse
import os
import sys
import time
from math import factorial
from collections import Counter

# ── Rotation group ────────────────────────────────────────────────────────────
# Face indices: U=0 L=1 F=2 R=3 B=4 D=5

# Three 90-degree generators (clockwise as seen from the indicated face):
X = (4, 1, 0, 3, 5, 2)  # R face: U->F->D->B->U  (L, R fixed)
Y = (0, 4, 1, 2, 3, 5)  # U face: F->R->B->L->F  (U, D fixed)
Z = (1, 5, 2, 0, 4, 3)  # F face: U->R->D->L->U  (F, B fixed)

def compose(a, b):
    return tuple(a[b[i]] for i in range(6))

# Generate all 24 rotation group elements via BFS
rots = [tuple(range(6))]
_seen = {tuple(range(6))}
_queue = [tuple(range(6))]
while _queue:
    _p = _queue.pop(0)
    for _g in (X, Y, Z):
        _np = compose(_g, _p)
        if _np not in _seen:
            _seen.add(_np)
            rots.append(_np)
            _queue.append(_np)

assert len(rots) == 24

# ── Key helpers ───────────────────────────────────────────────────────────────

def apply_rot(perm, key):
    return ''.join(key[perm[i]] for i in range(6))

def relabel(key):
    """Canonicalize digit labels by order of first appearance."""
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
    """Canonical form for reduction keys: lex-min over rotations with digit
    relabeling, restricted to rotations that preserve the canonical face-subset."""
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

def face_subset_of(red_key):
    return ''.join('0' if c != 'x' else 'x' for c in red_key)

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

# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dir', default='pipeline',
                   help='working directory for pipeline output '
                        '(default: pipeline, relative to this script)')
    p.add_argument('--max-cost', type=int, default=None, metavar='N',
                   help='filter bigseed.algs to cost <= N before expansion; '
                        'conjugations of surviving seeds may still exceed N '
                        'and are kept (not applied to costalgs output)')
    return p.parse_args()

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    root = os.path.dirname(os.path.abspath(__file__))

    # ── 1. Face-subset keys {0, x}^6 ─────────────────────────────────────────

    all_face_keys = set()
    for mask in range(64):
        key = ''.join('0' if (mask >> i) & 1 else 'x' for i in range(6))
        all_face_keys.add(canonical_subset(key))

    SKIP_FACE = {'xxxxxx', '0xxxxx'}
    face_keys = sorted(k for k in all_face_keys if k not in SKIP_FACE)

    print("Face-subset keys (for algutil -f; '0'=active, 'x'=fixed):")
    for k in face_keys:
        n = k.count('0')
        print(f"  {k}  ({n} active face{'s' if n != 1 else ''})")
    print()

    # ── 2. Reduction keys {0..5, x}^6 ────────────────────────────────────────

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

    filtered = sorted(
        (k for k in red_keys if 1 < space_size(k) <= 1_000_000_000),
        key=lambda k: (space_size(k), k)
    )

    # ── 3. Pipeline setup ─────────────────────────────────────────────────────

    os.chdir(root)
    if os.system("make") != 0:
        sys.exit("make failed")

    pipe_dir = args.dir if os.path.isabs(args.dir) else os.path.join(root, args.dir)
    os.makedirs(pipe_dir, exist_ok=True)
    os.chdir(pipe_dir)

    A = "../algutil"
    F = "../findpermnet"
    P = "../pareto"
    B = "../brobdicube"

    if args.max_cost is not None:
        print(f"Filtering bigseed.algs to max cost {args.max_cost:,} before expansion ...")
        print()

    # ── 3a. Expand seed algorithms ────────────────────────────────────────────

    print("Generating expanded.algs ...")
    if args.max_cost is not None:
        os.system(f"{A} -m {args.max_cost} < ../bigseed.algs 2>/dev/null"
                  f" | {A} -g -U 2>/dev/null | {A} -xu > expanded.algs 2>/dev/null")
    else:
        os.system(f"{A} -g -U < ../bigseed.algs 2>/dev/null | {A} -xu > expanded.algs 2>/dev/null")
    with open("expanded.algs") as f:
        total = sum(1 for _ in f)
    print(f"  {total:,} algorithms")
    print()

    # ── 3b. Generate costalgs for each face-subset key ────────────────────────

    print("Generating costalgs files for each face-subset key:")
    costalgs_counts = {}
    for fsk in face_keys:
        if fsk == '000000':
            os.system(f"{A} -c -d 1000 -pc < expanded.algs"
                      f" > expanded-{fsk}.costalgs 2>/dev/null")
        else:
            os.system(f"{A} -f {fsk} < expanded.algs > expanded-{fsk}.algs 2>/dev/null")
            os.system(f"{A} -c -d 1000 -pc < expanded-{fsk}.algs"
                      f" > expanded-{fsk}.costalgs 2>/dev/null")
        with open(f"expanded-{fsk}.costalgs") as f:
            n = sum(1 for _ in f)
        costalgs_counts[fsk] = n
        print(f"  expanded-{fsk}.costalgs: {n:,} algorithms")
    print()

    # ── 4. Run findpermnet for each reduction key ─────────────────────────────

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

    # ── 4b. Remove expanded* files (large; no longer needed) ─────────────────

    to_delete = ['expanded.algs']
    for fsk in face_keys:
        to_delete.append(f'expanded-{fsk}.costalgs')
        if fsk != '000000':
            to_delete.append(f'expanded-{fsk}.algs')
    deleted = 0
    for fname in to_delete:
        try:
            os.remove(fname)
            deleted += 1
        except FileNotFoundError:
            pass
    print(f"\nRemoved {deleted} expanded* intermediate file{'s' if deleted != 1 else ''}")

    # ── 5. Run pareto on all log files ────────────────────────────────────────

    algs_dir = "algs3"

    log_files = [f"permnet-{k}.log" for k in filtered]
    found_logs = [lf for lf in log_files if parse_permnet_log(lf)['found']]
    print(f"\nRunning pareto on {len(log_files)} log files"
          f" ({len(found_logs)} with solutions) -> {algs_dir}/")
    sys.stdout.flush()
    os.system(f"{P} -o {algs_dir} {' '.join(log_files)}")
    n_algs = sum(1 for f in os.listdir(algs_dir) if f.endswith('.alg'))
    print(f"  {algs_dir}/: {n_algs} .alg file{'s' if n_algs != 1 else ''}")

    # ── 6. Brobdicube planning sweep across a wide range of cube sizes ────────

    print(f"\nRunning brobdicube planning sweep: brobdicube -P 4 100000000000000000 -D {algs_dir}")
    sys.stdout.flush()
    t0 = time.time()
    os.system(f"{B} -P 4 100000000000000000 -D {algs_dir}")
    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s")

    os.chdir(root)


if __name__ == '__main__':
    main()
