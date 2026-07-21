#!/usr/bin/env python3
"""
Test pipeline: seed.algs -> expanded.algs / expanded.costalgs
               -> expanded-0x0x0x.algs / expanded-0x0x0x.costalgs
               -> permnet-phase1.log (all faces, key 010101)
               -> permnet-phase2.log (U/F/B only, key 0x1x2x)
               -> testalgs/ (Pareto-optimal .alg files)
               -> brobdicube -D testalgs -p 512 16384

Assumes seed.algs already exists in the parent directory (next to this script).
All intermediate files are written into a 'test/' subdirectory.
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))

os.chdir(root)
if os.system("make") != 0:
    sys.exit("make failed")

test_dir = os.path.join(root, "test")
os.makedirs(test_dir, exist_ok=True)
os.chdir(test_dir)

# Executables live one level up
A = "../algutil"
F = "../findpermnet"
P = "../pareto"
B = "../brobdicube"

# Step 1: greedy conjugation expansion -> plain algorithm list
os.system(f"{A} -g -U < ../seed.algs | {A} -xu > expanded.algs")

# Step 2: cost-annotate the full expanded set (include -pc for cycle notation)
os.system(f"{A} -c -d 1000 -pc < expanded.algs > expanded.costalgs")

# Step 3: filter to algorithms touching only U/F/B faces (key 0x0x0x)
os.system(f"{A} -f 0x0x0x < expanded.algs > expanded-0x0x0x.algs")

# Step 4: cost-annotate the filtered set (include -pc for cycle notation)
os.system(f"{A} -c -d 1000 -pc < expanded-0x0x0x.algs > expanded-0x0x0x.costalgs")

# Report line counts
for fname in ["expanded.algs", "expanded.costalgs",
              "expanded-0x0x0x.algs", "expanded-0x0x0x.costalgs"]:
    with open(fname) as f:
        print(f"{fname}: {sum(1 for _ in f)} lines")

# Step 5: phase 1 beam search — all-face algorithms, 2-color-group key
os.system(f"{F} -b 20 -w 20 -i 010101 < expanded.costalgs > permnet-phase1.log 2>&1")

# Step 6: phase 2 beam search — U/F/B-only algorithms, 3-face key
os.system(f"{F} -b 100 -w 100 -i 0x1x2x < expanded-0x0x0x.costalgs > permnet-phase2.log 2>&1")

# Step 7: collect Pareto-optimal solutions from both logs into testalgs/
os.system(f"{P} -o testalgs permnet-phase1.log permnet-phase2.log")

# Step 8: run brobdicube against the solution set
os.system(f"{B} -D testalgs -p 512 16384")

os.chdir(root)
