# actbeam6: Performance Optimization Notes

## Background

`actbeam5.cpp` is a beam search program that minimizes average element activations
for permutation network problems.  For the 12-zeros + 12-ones case (universe size
C(24,12) = 2,704,156 states), it was prohibitively slow even with beam width 20,
branching factor 20, and depth limit 63.

## Root Cause

The dominant cost was memory bandwidth from large array copies.  A `MaterializedState`
for the 24-element problem holds approximately:

| Field       | Size per state | Total (2.7M states) |
|-------------|----------------|---------------------|
| `best_mask` | 32 B (Mask256) | ~86 MB              |
| `best_cost` | 8 B            | ~22 MB              |
| `best_count`| 2 B            | ~5 MB               |
| `covered`   | 1/8 B          | ~338 KB             |
| `states`    | up to 8 B each | ~22 MB              |
| **Total**   |                | **~135 MB**         |

Two copy sites fired repeatedly:

1. **`extend_from_materialized` (line 1034 of actbeam5.cpp)**: did
   `MaterializedState out = st;` — a full 135 MB copy — for *every*
   (beam state, candidate) pair.  With W=20 and B=20, that is
   400 × 135 MB = **54 GB per depth step**, plus `apply_candidate_in_place`
   inside added another 108 MB snap copy per call (86 MB snap_mask +
   22 MB snap_cost), totalling 400 × 221 MB = **88 GB per step** just
   for the extend phase.

2. **`materialize` replay**: at depth d, the code called
   `apply_candidate_in_place` d times, each making a 108 MB snap.
   With W=20 beam states at depth 63, that is 20 × 63 × 108 MB ≈
   **136 GB per step** just in replay snap copies.

Grand total at the last depth step: roughly **230 GB of memcpy** per step,
repeated 63 times.

## Changes in actbeam6.cpp

### 1. `BeamState` — two new fields

```cpp
uint32_t parent_idx = 0;   // which beam_mats[] entry is this child's parent
int child_cand      = -1;  // candidate applied to produce this child
```

Used only during the depth-step transition; invisible to comparators and
dedup logic.

### 2. `ScoreScratch` struct (new)

Per-thread scratch arrays for `score_extension`:

```cpp
struct ScoreScratch {
    static constexpr uint16_t EMPTY = std::numeric_limits<uint16_t>::max();
    std::vector<uint16_t> temp_count;  // EMPTY = not touched this call
    std::vector<uint64_t> temp_cost;
    std::vector<uint32_t> touched;     // ranks modified; reset before return
    void init(uint32_t total);
};
```

Allocated once per thread per depth step; reused across all candidate calls.

### 3. `score_extension` function (new, ~80 lines)

Replaces `extend_from_materialized`.  Reads the parent `MaterializedState`
as `const`; tracks virtual destination-rank updates in `ScoreScratch`;
resets scratch before returning.  Returns a `BeamState` with scores and
`parent_idx`/`child_cand` set.

Key correctness point: because sources are always read from the immutable
parent (never from a modified copy), the snapshot that `apply_candidate_in_place`
requires to prevent false transitive chains is not needed here.

Approximate tie-break: `b.tie = st.tie ^ (cand_idx * 0x9e3779b97f4a7c15ULL)`.
The exact tie is recomputed during survivor materialization (phase 2) and
corrected before the next depth step.

### 4. `beam_mats` — cached MaterializedStates (main loop change)

A `std::vector<MaterializedState> beam_mats` is maintained in parallel with
`beam`.  `beam_mats[i]` is the fully computed DP state for `beam[i]`.

The main loop now has two phases per depth step:

**Phase 1 — scoring (parallel, read-only, no large copies)**

Each thread processes its assigned slice of `beam_mats[begin..end)`.  For
each beam state, for each candidate, `score_extension` is called.  No
MaterializedState is copied; no snap arrays are allocated.

**Phase 2 — survivor materialization (sequential in main thread)**

After global merge and trim, for each surviving child:

```cpp
MaterializedState mat = beam_mats[b.parent_idx];  // one copy per survivor
apply_candidate_in_place(mat, b.child_cand, depth, ...);
b.tie = mat.tie;   // correct the approximate tie
next_mats.push_back(std::move(mat));
```

Cost: W copies of ~135 MB, down from W×B×221 MB.
For W=B=20 this is a 20× reduction in copy volume for the extend phase,
plus the full replay (materialize-from-scratch) is eliminated.

### 5. `extend_from_materialized` removed

Dead code; deleted.

## Memory cost of beam_mats

Storing W `MaterializedState` objects simultaneously:

- W=20, n=24: 20 × 135 MB = **2.7 GB**
- During phase 2 (before old mats are freed): peak ≈ 4–5 GB

This is the main memory trade-off of this approach.

## Measured Speedup

### Test utility files

| File | Description |
|------|-------------|
| `/tmp/perms8.txt`  | All 40,320 permutations of {0..7}, image form, generated with `itertools.permutations` |
| `/tmp/perms12.txt` | 10,010 permutations of {0..11}: 11 adjacent transpositions + 10,000 random distinct |
| `/tmp/perms16.txt` | 5,000 random distinct permutations of {0..15} |
| `/tmp/perms20.txt` | 3,000 random distinct permutations of {0..19} |

All generated with Python one-liners (see session history); not committed.

### Results

Comparing `actbeam5` vs `actbeam6`, seed 42, `-w 20 -b 20`:

| n  | Universe   | Mat size | Depth | actbeam5 | actbeam6 | Speedup |
|----|-----------|----------|-------|----------|----------|---------|
| 8  | 70        | ~2 KB    | 20    | 0.28 s   | 0.25 s   | 1.1×    |
| 12 | 924       | ~38 KB   | 40    | 5.42 s   | 4.40 s   | 1.2×    |
| 16 | 12,870    | ~542 KB  | 10    | 1.30 s   | 1.17 s   | 1.1×    |
| 20 | 184,756   | ~11 MB   | 6     | 2.46 s   | 0.80 s   | **3.1×** |

For n ≤ 16 the MaterializedState fits in L3 cache and the speedup is minimal.
At n=20 (~11 MB) the data exceeds L3 and the 3× speedup is clearly visible.
The n=24 case (135 MB, depth 63) should yield a far larger speedup (estimated
20–50×) since both the state-space size (14.6×) and depth (10.5×) scale the
memory pressure far beyond what was tested here.

## actbeam7: Intermediate Solutions and Early Exit

### Changes

**1. Intermediate solution printing**

Previously the program printed nothing until the depth limit was reached, then
printed at most one solution.  `actbeam7` prints a `FOUND` block at every depth
where the beam first achieves full coverage (`cover_count == ranker.total`),
using the already-cached `beam_mats[0]` — no extra replay or materialization
cost.

**2. Early exit when activations plateau**

After printing each complete solution the program compares `sum_activations`
against the previous complete solution:

```cpp
if (have_solution && best.sum_activations >= prev_sum_activations) {
    return 0;   // no improvement — stop
}
prev_sum_activations = best.sum_activations;
have_solution = true;
```

The exit fires when activations are equal (or, in pathological beam-search
cases, worse).  This prints one "confirming" solution at the plateau depth
then exits, rather than running indefinitely.

**3. `print_solution` refactored**

Split into a core overload taking `const MaterializedState&` (used by the
main loop with the cached mat) and the materialize-then-print wrapper retained
for any future callers.  The `materialize` function itself became dead code
and was removed.

### Verification

Test input `/tmp/swaps12.txt`: all C(12,2) = 66 single-swap permutations of
{0..11} (generated with a Python two-loop one-liner; not committed).  These
have limited coverage ability and cause the beam to plateau quickly.

```
$ ./actbeam7 -w 20 -b 20 -e 40 -i 000000111111 -s 42 < /tmp/swaps12.txt
```

Selected output (stderr suppressed):

```
FOUND depth 17  sum_activations 3653
FOUND depth 18  sum_activations 3450
...
FOUND depth 37  sum_activations 2773
FOUND depth 38  sum_activations 2773   ← equal to depth 37, exit triggered
```

The program exited after depth 38 rather than running to `-e 40`.

## actbeam8: Minimize Max Activations, then Average

### Objective change

The primary optimization criterion shifted from minimizing the *average* number
of activations to minimizing the *worst case* (maximum over all target states),
with average as a secondary tiebreaker.  New beam sort order:

```
cover_count ↓  →  max_activations ↑  →  sum_activations ↑  →  sum_activation_cost ↑
```

### Changes

**`MaterializedState` and `BeamState`**: added `uint16_t max_activations`.

**`make_start_materialized`**: initializes `max_activations = 0` (the single
starting state has an empty activation mask, so its count is 0).

**`score_extension` (upper-bound approximation)**:
Exact max tracking in the read-only scoring path would require an O(ranker.total)
scan per call — at W×B = 400 calls per depth step on the 2.7M-state problem
that would dominate runtime.  Instead, `max_activations` is only updated upward:
newly covered states can raise it, and improvements to already-covered states
also raise it if the new count is higher (which cannot happen by construction,
but the code is written conservatively).  Crucially, the max is never decreased
in `score_extension`, so it is a guaranteed upper bound on the true child max.
The only impact on quality is that a beam state whose max genuinely decreased
(because the unique worst-case state found a shorter path) looks the same as one
whose max did not change; the correct value is restored in phase 2.

**`apply_candidate_in_place` (exact max)**:
After the main DP loop, one O(covered_count) scan recomputes the exact maximum:

```cpp
out.max_activations = 0;
for (uint32_t r = 0; r < ranker.total; ++r)
    if (test_bit(out.covered, r) && out.best_count[r] > out.max_activations)
        out.max_activations = out.best_count[r];
```

Called only W times per depth step (phase 2 survivor materialization), so
the scan cost is acceptable.

**Phase 2 survivor loop**: after `apply_candidate_in_place`, both `b.tie`
and `b.max_activations` are copied from the exact materialized result,
correcting the upper-bound estimate carried through from `score_extension`.

**`beam_dedup_key`**: updated to include `max_activations`.

**`print_solution`**: now prints `max_activations` immediately after `coverage`.

**Progress line**: now shows `maxact=` alongside `sumact=` and `avgact=`.

**Early-exit condition**: extended to require that *both* `max_activations`
and `sum_activations` fail to improve before exiting:

```cpp
if (have_solution &&
    best.max_activations >= prev_max_activations &&
    best.sum_activations >= prev_sum_activations) {
    return 0;
}
```

This lets the search continue if only the average is still improving (max
plateau) or only the max is still improving (unlikely but possible).

### Verification

Using `/tmp/swaps12.txt` (66 single-swap permutations of {0..11}, `-w 20 -b 20
-e 80 -i 000000111111 -s 42`):

```
FOUND depth 17  max_activations 7   sum_activations 4015
FOUND depth 18  max_activations 7   sum_activations 3713
FOUND depth 19  max_activations 6   sum_activations 3625
...
FOUND depth 43  max_activations 6   sum_activations 2772
FOUND depth 44  max_activations 6   sum_activations 2772   ← both plateau, exit
```

max stabilised at 6 from depth 19 onward; sum plateaued at 2772 at depth 44,
triggering the combined exit condition.

## addmax.cpp: Back-fill max_activations into actbeam7 Logs

### Purpose

Existing actbeam7 log files contain complete `FOUND` blocks but lack the
`max_activations` line that actbeam8 adds.  `addmax.cpp` is a standalone
post-processor that reads an actbeam7 stdout log from stdin, computes the
missing `max_activations` for each block, and writes the modified log to
stdout — with the new line inserted after the `coverage N/T` line, exactly
as actbeam8 would have printed it.

The program needs neither the original permutation input file nor any
command-line arguments; all information required for the DP replay is already
present in each FOUND block (`initial XXXXX` and the per-candidate cycle
notations).

### Key design decisions

**Simplified DP (no Mask256)**

The full DP in actbeam8 tracks a 256-bit activation mask per rank to support
exact cost tie-breaking.  For computing `max_activations` only the hop count
matters: the minimum number of elements (sequence positions) needed to reach
each rank.  Because each element at sequence position `i` sets a unique bit
`i` that was never set before, `new_count = snap_count[source] + 1` always
holds.  So the per-rank state reduces from `Mask256 + uint16_t + uint64_t`
to a single `uint16_t best_count[r]`, and the Mask256 layer is omitted
entirely.

**Snapshot fix retained**

The snapshot of source counts before each element's loop is still required
to prevent false transitive chains: if destination rank `r` is also a source
state later in the same loop, an in-place update to `best_count[r]` would
corrupt subsequent iterations.  The snapshot cost is `O(ranker.total)` per
element, which is acceptable for a post-processing tool.

**Codec/Ranker caching**

The expensive `Ranker::init` call (which builds the full rank tables) is
cached across consecutive FOUND blocks with the same initial string.  For
logs that all share one problem instance (the common case) it fires only once.

**Pass-through for actbeam8 logs**

If the line after `coverage N/T` already begins with `max_activations `, the
block is emitted unchanged.  This makes `addmax` safe to run on logs that are
already in actbeam8 format or on mixed logs.

### Changes relative to actbeam8 infrastructure

| actbeam8 component         | addmax treatment             |
|---------------------------|------------------------------|
| `Mask256`, mask ops        | omitted                      |
| `BeamState`, `MaterializedState` | omitted             |
| `ScoreScratch`             | omitted (no parallel scoring)|
| `ThreadLocalLUTCache`      | omitted (LUTs built per step)|
| `Candidate`, `CostBucket`  | omitted                      |
| `Codec`, `ChunkLUT`, `Ranker` | copied verbatim           |
| `build_chunk_lut`, `apply_lut` | copied verbatim          |
| `set_bit`, `test_bit`, `symbol_value` | copied verbatim  |

### Verification

Generated a reference actbeam8 log, stripped its `max_activations` lines to
simulate an actbeam7 log, ran `addmax`, and diffed against the original:

```
$ ./actbeam8 -w 10 -b 10 -e 25 -i 000000111111 -s 42 < /tmp/swaps12.txt 2>/dev/null \
    > /tmp/ab8_log.txt
$ grep -v '^max_activations ' /tmp/ab8_log.txt > /tmp/ab8_stripped.txt
$ ./addmax < /tmp/ab8_stripped.txt > /tmp/ab8_restored.txt
$ diff /tmp/ab8_log.txt /tmp/ab8_restored.txt
(no output — exact match)
```

Also verified that running `addmax` on an already-processed actbeam8 log
produces bit-for-bit identical output (pass-through path).

## actbeam9.cpp: Sum-of-Squares Activation Scoring

### Objective change

Instead of minimizing max activations then average, actbeam9 minimizes the
**sum of squares** of per-state activation counts as the primary criterion,
with sum (average) as a secondary tiebreaker.  The intuition is that squaring
penalises outliers more heavily than summing, naturally driving the worst-case
count down without an explicit max term.

New beam sort order:
```
cover_count ↓  →  sum_sq_activations ↑  →  sum_activations ↑  →  sum_activation_cost ↑
```

### Changes

- `uint64_t sum_sq_activations` added to both `BeamState` and `MaterializedState`.
- `score_extension`: tracks `sum_sq_activations` alongside `sum_activations`
  (newly covered: `+= count²`; improvement: `−old² + new²`).
- `apply_candidate_in_place`: same incremental sq tracking; exact max rescan
  retained at end.
- `beam_better`: replaced `max_activations` primary with `sum_sq_activations`.
- `beam_dedup_key`: `sum_sq_activations` replaces `max_activations` in the hash.
- Phase 2 survivor loop: corrects `b.sum_sq_activations` from the exact mat.
- `print_solution`: emits `sum_sq_activations` after `max_activations`.
- Progress line: adds `sumsq=` field.
- Early exit: fires when `sum_sq_activations` fails to improve (single
  criterion, replacing the dual max+sum condition from actbeam8).
- `max_activations` is still tracked and reported but no longer drives search.

## actbeam10.cpp: Eliminate best_mask from Hot DP Path

### Motivation

For n=24 (C(24,12)=2,704,156 states), each `MaterializedState` held a
`std::vector<Mask256> best_mask` of 2.7M × 32 bytes = **86 MB**.  This
dominated memory and memory-bandwidth costs:

| Site | Old cost | New cost |
|------|----------|----------|
| `beam_mats` storage (W=20) | 20 × 135 MB = 2.7 GB | 20 × 49 MB = 980 MB |
| `snap_mask` copy per phase-2 call | 86 MB | 0 |
| `snap_count` copy per phase-2 call | 0 | 5.4 MB |
| `score_extension` mask read per source state | 32 B (Mask256) | 2 B (uint16_t) |

### Key insight

`new_count = popcount(best_mask[r0] | elem_bit) = best_count[r0] + 1` always,
because each element at sequence position `i` sets bit `i`, a bit that has
never appeared in any mask before position `i`.  So `best_mask` is never
actually needed to compute counts or costs — only to reconstruct which elements
activated which states when printing.

### Changes

**`MaterializedState`**: `best_mask` field removed entirely.

**`score_extension`**: `elem_bit` and all `mask256_*` calls removed; replaced
with `new_count = st.best_count[r0] + 1`.  `elem_idx` parameter dropped (now
unused).

**`apply_candidate_in_place`**: `elem_bit`, `snap_mask` (86 MB copy), and all
`best_mask[r]` writes removed.  A `snap_count` copy (5.4 MB) replaces
`snap_mask` to prevent false transitive chains on counts.  `elem_idx` parameter
dropped.

**`compute_activation_counts`** (new): called only when printing, replays the
winning sequence from scratch with full mask tracking to recover per-element
activation counts.  Uses a partial snapshot (only source-rank entries, not the
full 86 MB array) so snapshot cost is O(states.size()) per step.  Allocates
86 MB temporarily for the replay; freed on return.

**`print_solution`**: accepts `ThreadLocalLUTCache&`; calls
`compute_activation_counts(st.seq, ..., mat.states[0])` to obtain
`activation_counts[]` before printing.

**Progress line**: `elapsed=Xs` field added showing wall-clock time for each
depth step (time since the previous depth line, not total elapsed time).

### Memory reduction summary

| Component | actbeam9 | actbeam10 |
|-----------|----------|-----------|
| `best_mask` per state | 86 MB | 0 |
| `best_count` per state | 5.4 MB | 5.4 MB |
| `best_cost` per state | 21 MB | 21 MB |
| `covered` per state | 0.34 MB | 0.34 MB |
| `states` per state | up to 21 MB | up to 21 MB |
| **MaterializedState total** | **~135 MB** | **~49 MB** |
| **beam_mats (W=20)** | **~2.7 GB** | **~980 MB** |

## actbeam10.cpp: Operational Improvements

### argv logging

The first line of stdout is now `- argv[0] argv[1] ...` — a hyphen followed by
the full argument list.  This allows log files to self-document how they were
produced.

### Per-phase timing

The `elapsed=Xs` field is replaced with `p1=Xs p2=Ys`, where:
- `p1` = wall time for Phase 1 (scoring all candidate branches)
- `p2` = wall time for Phase 2 (materializing survivors)

### Phase 2 parallelization

Survivor materializations are independent (each child reads from its parent's
`beam_mats` entry and writes to its own `next_mats` slot), so phase 2 is now
parallelized across the same `-t N` worker threads as phase 1:

```cpp
std::vector<MaterializedState> next_mats(next.size());
auto p2_ranges = make_ranges(next.size(), num_threads);
// each thread gets its own ThreadLocalLUTCache and writes to distinct indices
for (auto& th : p2_workers) th.join();
```

This addresses the "Phase 2 parallelization" item from the previous remaining
opportunities list.

### Async print_solution

`compute_activation_counts` — called every depth step once fully covered —
does a full Mask256 replay over all 2.7M states, allocating ~172 MB
(`best_mask` + `snap_mask`).  For n=24 this takes several seconds and was the
dominant single-threaded bottleneck.

Each found solution now launches a background `std::thread` that runs
`print_solution` independently.  The depth loop continues immediately.

Key design points:
- `print_solution` now takes `BeamState` by value (copied at thread-launch
  time), `State init_state` by value, `size_t lut_budget`, and
  `std::mutex& output_mutex`.
- All scalar stats (`cover_count`, `max_activations`, `sum_sq_activations`,
  `sum_activations`, `sum_activation_cost`) come from the `BeamState` copy,
  which is already exact after phase 2.  No `MaterializedState` is needed.
- `init_state` (`beam_mats[0].states[0]`) is constant across all depths and
  captured once before the depth loop.
- `output_mutex` serialises stdout if two print threads happen to finish
  concurrently.
- `print_threads` is joined at every exit point (three in actbeam10, four in
  metalbeam1) to ensure clean shutdown.

### `-q` quick print option

With `-q`, found solutions are printed immediately from the main thread with
no background thread and no `compute_activation_counts` call.  Output contains
all scalar stats from `BeamState` plus each step's candidate id, cost, and
permutation in cycle notation — but no per-step activation counts or
`contribution_num` fields.  This is lossless for the information already
available in `BeamState`; the per-step activation counts can be reconstructed
offline by a separate tool.

## metalbeam1.cpp: GPU-Accelerated Scoring via Apple Silicon Metal

### Motivation

Phase 1 scoring iterates over all covered states for every (beam-state, candidate)
pair.  With beam width 256 and branching factor 256, this is 65,536 candidates ×
2.7M states = 178 billion operations per depth step on a single-threaded CPU.
Apple Silicon's GPU (unified memory, no transfer overhead) can score all
candidates in parallel with one dispatch.

### Files created

| File | Purpose |
|------|---------|
| `beam_score.metal` | MSL compute shader (also embedded in metal_score.mm) |
| `metal_score.h` | Pure C++ interface; no ObjC headers |
| `metal_score.mm` | ObjC++ implementation of MetalScorer |
| `metalbeam1.cpp` | Fork of actbeam10.cpp with Metal Phase 1 |

Build command:
```
clang++ -O3 -g -std=c++17 metalbeam1.cpp metal_score.mm \
    -framework Metal -framework Foundation -o metalbeam1
```

### GPU shader design (`beam_score.metal`)

**One threadgroup per candidate; 512 threads per threadgroup.**  Each thread
strides over states `[tid, tid+512, tid+1024, ...]`, accumulates deltas
(`delta_sum_sq`, `delta_sum_act`, `delta_sum_cost`, `delta_cover`) in
registers, then reduces into a `threadgroup long[512]` array.  Thread 0
writes the final `ScoreOut` struct for that candidate.  No atomics; no write
conflicts.

**ChunkLUT on GPU**: the n=24 LUT decomposes into 2 chunks of 4096 entries
each (64 KB per candidate), stored in a flat `device const ulong* lut_data`
buffer.  Candidate `i`'s entries start at `lut_offsets[i]`.  With 256
candidates that is 16 MB — well within GPU device memory.

**Ranker on GPU**: the two-level rank tables (`high_base`, `high_low_table`,
`low_rank_flat`) are copied once at startup.  For n=24: 4096 + 4096 + 13×4096
= ~80 KB total.  `gpu_rank` mirrors the C++ `Ranker::rank()` exactly.

**Approximation**: the GPU ignores within-permutation double-counting (two
covered sources mapping to the same destination).  This is a valid
approximation for ranking candidates; the winner is materialised exactly on CPU
in phase 2.  `out_max_act` is set conservatively to `base_max_act`; the true
max is recomputed by `apply_candidate_in_place`.

### metal_score.h / metal_score.mm interface

```cpp
// Pure C++ — no ObjC headers leak into metalbeam1.cpp
MetalScorer* metal_scorer_create(uint32_t n_states,
    const uint64_t* state_of_rank,
    const uint32_t* high_base, const uint32_t* high_low_table, uint64_t high_mod,
    const uint32_t* low_rank_flat, uint64_t low_mod, uint32_t n_low_tables,
    uint64_t split_mod, const uint32_t* chunk_mods, int nchunks);

void metal_scorer_score(MetalScorer* ms,
    const uint16_t* best_count, const uint64_t* best_cost,
    const uint64_t* covered,
    uint64_t base_sum_sq, uint64_t base_sum_act, uint64_t base_sum_cost,
    uint32_t base_cover_count, uint16_t base_max_act,
    const uint64_t* lut_flat, const uint32_t* lut_offsets,
    const uint64_t* cand_costs, int n_candidates,
    uint64_t* out_sum_sq, uint64_t* out_sum_act, uint64_t* out_sum_cost,
    uint32_t* out_cover_count, uint16_t* out_max_act);

void metal_scorer_destroy(MetalScorer* ms);
```

`MetalScorer` holds persistent `id<MTLBuffer>` objects for read-only data
(state_of_rank, rank tables, chunk_mods) that are uploaded once at creation.
Per-call buffers (best_count, best_cost, covered, lut_flat, cand_params,
out_scores) are grown on demand via `ensure_buf`.  Apple Silicon unified memory
means `memcpy` into the buffer is the only "transfer" cost.

### metalbeam1.cpp changes vs actbeam10.cpp

- `#include "metal_score.h"` at top; `metal_scorer_destroy` at all exit paths.
- After `ranker.init()`: flatten `low_rank_tables` into `low_rank_flat` and
  call `metal_scorer_create`.
- Phase 1 replaced with GPU path (when Metal is available): build flat
  `lut_flat` array (n_candidates × lut_entries_per_cand uint64 entries),
  call `metal_scorer_score`, construct `BeamState` children from output arrays.
  CPU multi-threaded fallback retained.
- Phase 2 parallelised identically to actbeam10.
- All operational improvements from actbeam10 (argv print, p1/p2 timing, async
  print, `-q` flag) are present.

### Expected speedup

For `-b 256` (256 candidate branches per beam state), the GPU dispatches 256
threadgroups simultaneously, each scoring 2.7M states in ~0.5–2 ms (M-series
GPU has ~100–400 GB/s memory bandwidth).  CPU single-threaded Phase 1 at the
same scale would take ~20–200 s.  Practical speedup depends on how many
candidates fit in the trial set and on beam width; the GPU path is most
effective with large `-b` values that are impractical on CPU.

### Wall-clock elapsed timer and read timer

`program_start` is captured at the very first line of `main()` (before argument
parsing or input reading).  Each depth progress line now includes:

- `elapsed=Xs` — seconds since `program_start` to the end of that depth step,
  giving an immediate diagnostic for startup/read overhead at depth 1.
- `read=Xs` — seconds spent in `read_perms()` specifically, printed once in
  the startup banner line.

Both fields are in actbeam10 and metalbeam1.

### C-style I/O rewrite for read_perms

`std::getline` / `std::string` parsing was replaced with `fgets` +
`const char*` throughout `read_perms`, `parse_perm_line`,
`leading_cost_value`, and `parse_ints`.  Highlights:

- `parse_ints(const char* p, const char* end)` — scans with `strtol`,
  no `substr` or `stoi`.
- `leading_cost_value(const char* line)` — pointer walk + `strtoull`
  with `endptr`; no string construction.
- `parse_perm_line` — pure pointer arithmetic; `buf.data()` passed directly.
- Lines >65535 bytes detected via the fgets sentinel byte, not `line.back()`.
- `std::string text` removed from `Candidate` struct (was unused after
  `parse_perm_line` was rewritten).

### stderr → stdout migration

All progress/informational output (`depth N:` lines, `NOT FOUND`, startup
banner, etc.) migrated from `std::cerr` to `std::cout`.  Only genuine
exceptions caught in the outer `try/catch` still go to `std::cerr`.  This
lets a single stdout redirect capture the complete log.

## collect_solutions.cpp: Pareto-Optimal Solution Collector

`collect_solutions.cpp` is a standalone post-processor that collects FOUND
blocks from any number of actbeam10/metalbeam1 log files, computes a
Pareto front, and writes the winning solutions as `.alg` files.

### Usage

```
./collect_solutions [-A ALGSDIR] [-o OUTDIR] [-d MAXDEPTH] [-n] LOGFILE...
```

- `-A ALGSDIR` — directory containing `bycosth0cx-{dotkey}.algs` move files
  (default: `/Users/rokicki/twsearch2`).  Move strings are loaded on demand,
  one file per dot-key, and cached for the run.
- `-o OUTDIR` — output directory (default: `algs3`).
- `-d MAXDEPTH` — discard solutions deeper than MAXDEPTH.
- `-n` — dry run: count how many solutions would be written without
  creating any files or the output directory.
- `LOGFILE...` — one or more log files; **stdin is not supported** (log
  filenames must be supplied explicitly so the dot-key can be extracted).

### Dot-key and shortkey

Each log filename embeds a 6-character **dot-key** (e.g. `.....0`, `0.0.00`)
identifying which faces were active during that beam search run.  The tool
extracts this key and combines it with the solution's `initial` string to
derive a 6-character **brob shortkey** (e.g. `01111x`, `x0x1xx`) — the face
partition identifier used by brob3.

Rules:
- A `.` in the dot-key means that face was active.  The four elements of that
  face in `initial` are all the same digit; that digit becomes the shortkey
  character for that face.
- A digit (`0`/`1`/`2`) in the dot-key means the face was already solved and
  fixed; that position becomes `x` in the shortkey.

The special all-active dot-key `......` loads `bycosth0cx.algs`; all other
dot-keys load `bycosth0cx-{dotkey}.algs`.

### Parsing

State-machine parser handles both full and `-q` format FOUND blocks.
Fields extracted per block:

| Field | Source line |
|-------|-------------|
| depth | `FOUND depth N` |
| coverage | `coverage M/T` |
| max_activations | `max_activations N` |
| sum_activations | `sum_activations N` |
| average_activations | `average_activations X` |
| initial | `initial XXXX` |
| candidates + cost fields | `N: candidate M cost COST ...` lines |

Only full-coverage blocks (`M == T`) are kept.  Each candidate's `cost`
token is decomposed into three independent fields of three decimal digits each
(`cost/1000000`, `(cost/1000)%1000`, `cost%1000`), accumulated separately in
`sum_cost_f1/f2/f3` to prevent field overflow at large depths.

### Algorithm

1. Parse all log files, collecting full-coverage solutions.
2. Compute the brob shortkey for each solution from the log filename.
   Drop any solution whose shortkey cannot be determined.
3. Deduplicate by **(shortkey, candidate sequence)** — candidate indices are
   pool-specific, so the same index sequence from different dot-keys is a
   different algorithm.
4. Compute Pareto front **per shortkey** on three metrics (all less-is-better):
   depth, max_activations, sum_activations.  Solutions from different shortkeys
   represent different permutation pools and are never compared.
5. Sort by (shortkey, depth, sum_activations, max_activations).
6. For each non-dominated solution:
   - Load the appropriate `.algs` move file for this dot-key (cached).
   - Build alg content: one space-prefixed move string per candidate line.
   - Write to `OUTDIR/{shortkey}-{depth:04d}-{f1:08d}-{f2:08d}-{f3:08d}-{hash}.alg`.
   - Print a `SOLUTION ... shortkey=... initial=...` summary to stdout.
7. Statistics go to stderr.

### Notes

- The FNV-1a hash of the candidate sequence provides a stable unique filename
  suffix even after deduplication.
- Multiple dot-keys in one run are fully supported; each produces its own
  Pareto front keyed by shortkey.

## brob3.cpp: Algorithm Application with algs3 Support

`brob3.cpp` is derived from `brob2.cpp` with the following changes:

### algs3 directory and dual format support

`readalgs()` reads from `./algs3` instead of `./algs2`.  It accepts files in
two formats, distinguished by counting hyphens in the filename:

| Hyphens | Format | Depth filter |
|---------|--------|-------------|
| 2 | old: `{key}-{packed_cost}-{hash}.alg` | none |
| 5 | new: `{key}-{depth:04d}-{f1}-{f2}-{f3}-{hash}.alg` | depth > `planmask::max_stages()` skipped |

Any other hyphen count is silently skipped.

### `planmask` type

The `int` plan bitmask that flowed between `plan()`, `getplan()`, `doplan()`,
and the `movemasks` array has been replaced with an explicit `planmask` struct:

```cpp
struct planmask {
   uint64_t w = 0 ;
   static planmask none() ;          // sentinel (all bits set)
   static int max_stages() ;         // 63: capacity of current representation
   bool is_none() const ;
   planmask set(int i) const ;       // return copy with bit i set
   bool test(int i) const ;          // test bit i
   int popcount() const ;
   bool operator==(planmask) const ;
   bool operator!=(planmask) const ;
} ;
```

This removes the N=31 limit (int had 31 usable bits; uint64_t has 63 usable
bits with the all-ones sentinel reserved).  The interface is designed so that
when STAGE > 63 is needed, only the struct internals change — all call sites
use `set`/`test`/`popcount`/`is_none`.

`max_stages()` is used directly in `readalgs()` to set the depth filter,
so it automatically tracks any future capacity change.

### Mx ≥ My invariant and assert

`main()` previously swapped the Mx/My assignment relative to `dosolve()`,
causing `phasewidth > cube.My` on odd-N inputs.  This triggered a latent bug
in `wplan::plan()` where stale count values from a prior `My` were not reset
before re-optimising for the capped (smaller) `My`, resulting in a heap
buffer overread.

Fixed by aligning `main()`'s computation with `dosolve()`:

```cpp
ll My = (N - 2) / 2 ;   // smaller half
ll Mx = (N - 2) - My ;   // larger half (Mx >= My always)
```

This same fix was applied to `brob2.cpp`.

An assert was added at the end of `wplan::plan()` to catch any future
recurrence:

```cpp
assert(count1 * w1 + count2 * w2 == My) ;
```

### Dual-mode planning: greedy and optimal

Each `.alg` file loaded by `readalgs()` now creates **two** `netgen` objects —
one greedy, one optimal — which are both added to `algset`.  `bestcost()` tries
all candidates (as before) and naturally picks whichever mode gives the lower
overall cost.  `explain()` labels each chosen step `[greedy]` or `[optimal]`.

**Greedy mode** (current algorithm, unchanged): one pass through the algorithms
in index order.  For each algorithm `i`, every already-planned state `j` is
extended; the new state `d = alg_i(j)` is recorded if no plan exists for it
yet, or if the new plan has fewer algorithms (smaller popcount).  Plans are
committed at the end of each round and never revised.  Result: each state gets
the plan with the fewest algorithms.

**Optimal mode**: same BFS structure, two differences:

1. Comparison uses total algorithm cost `sum(algs[i].cost for i in plan)`
   instead of popcount.
2. The commit step updates a state's plan even if one already exists, whenever
   the new plan is cheaper.

Because both modes commit at the **end** of each round (not during the inner
loop), `plans[j].test(i)` is never true when processing round `i` — the same
invariant that makes the greedy plan correct also applies here, with no risk of
recording the same algorithm bit twice.  The optimal mode is a single-pass
cost-greedy heuristic; it is not provably globally optimal, but it produces
different `stageprob` distributions than greedy and generally prefers plans
with fewer but cheaper algorithms.

Prob cache files use mode-specific suffixes so the two sets never collide:
`{filename}.prob-greedy` and `{filename}.prob-optimal`.

When a rotated copy of a `netgen` is needed for `explain()`, the `optimal` flag
is propagated to the copy so it plans in the correct mode.

Build command (unchanged):
```
g++ -march=native -Wall -Wsign-compare -g -std=gnu++20 -O3 -o brob3 brob3.cpp
```

### New options

**`-D dirname`** — directory to read algorithm files from and write prob caches
to (default: `./algs3`).  Previously hardcoded; now a global `algsdir` set by
the flag.

**`-P start end`** — planning-only sweep mode.  Runs `bestcost` + `explain` for
sizes `start, start*2, start*4, ...` up to `end`, then exits.  No solving is
done.  Other flags (`-w`, `-p`, etc.) compose with `-P` as usual.  Caches
(`bestcostcache`, `bestcostplan`) are cleared between sizes so each N is
planned independently.  Example: `./brob3 -P 5 100 -w` plans worst-case costs
for N = 5, 10, 20, 40, 80.

## brob4.cpp: Extended planmask with Compact Mode

`brob4.cpp` is derived from `brob3.cpp` with a single structural change: the
`planmask` type gains a dual-mode representation so that algorithm sets with
more than 63 stages can be handled.

### Dual-mode planmask

| Mode | Trigger | Representation |
|------|---------|----------------|
| Bitmask | MSB (bit 63) = 0 | Bits 0..62 = which algorithms to apply |
| Compact | MSB (bit 63) = 1 | 7 sorted 9-bit fields; value 511 = empty slot |

**Bitmask mode** is identical to brob3's planmask.  Algorithm indices 0..62
map directly to bits.  `max_stages()` previously returned 63; now it returns
511 (the highest index storable in compact mode).

**Compact mode** stores up to 7 algorithm indices (0..510) in 7 × 9 = 63 bits,
always sorted ascending, with 511 filling unused slots.  The MSB set and all
fields = 511 equals `~0ULL` — the "not reached" sentinel — but this is not a
valid plan since compact mode is only entered by adding at least one algorithm,
so the two are never confused.

**Transition**: `set(i)` detects the mode:
- Bitmask + i < 63 → set bit i (stay bitmask).
- Bitmask + i ≥ 63 → extract set bits as a sorted compact sequence, then
  `compact_insert(i)`.  If popcount ≥ 7 (would need 8 slots), returns `none()`
  to signal overflow.
- Compact + any i → `compact_insert(i)`.  Returns `none()` if all 7 slots are
  full.

**Overflow handling**: both `plan_greedy()` and `plan_optimal()` check that
`set()` did not return `none()` before updating `next[d]`.  This silently skips
paths that would require more than 7 compact-mode activations, retaining any
existing (bitmask) plan for the state.  The user's expectation is that the first
63 algorithms cover all states, so compact mode is used only for cost
improvement, not initial coverage.

**`doplan()` change**: branches on `is_bitmask()`.  Bitmask path unchanged
(iterate bits 62..0).  Compact path iterates fields 6..0 descending,
skipping 511, and applies `algs[pi].perm` for each non-empty field — equivalent
order since fields are sorted ascending (highest index is in field 6).

Build:
```
g++ -march=native -Wall -Wsign-compare -g -std=gnu++20 -O3 -o brob4 brob4.cpp
```

### Startup speed: MINWRITE and prob2 cache invalidation

Two issues caused slow subsequent startups even after `.prob` files were written:

**1. MINWRITE threshold lowered from 65536 to 8192.**  `ensureprobs()` only writes a
`.prob2` cache file when `goal > MINWRITE`.  At 65536, many medium-goal netgens
(e.g. goal = 34650) were always re-planned from scratch on every run.  Lowering
to 8192 causes those netgens to be cached on the first run and skipped on all
subsequent runs.

**2. Cache suffix renamed to `.prob2-greedy` / `.prob2-optimal`** (from `.prob-greedy`
/ `.prob-optimal`).  The rename invalidates stale caches written with an earlier
version of the planner; any file with the old suffix is silently ignored.

**3. `f->plan()` removed from `explain()`.**  `explain()` was calling `plan()` on
the netgen even when a warm prob cache existed.  Since `getplan()` already
re-plans lazily when `plans` is empty, the explicit call was redundant.

### Cost encoding: base-10000 fields

`blockalgo::cost` is a single `int64_t` encoding four fields `(v, h, o, c)`:

```cpp
cost = c + 10000 * (o + 10000 * (h + 10000LL * v)) ;
```

Each field occupies four decimal digits, allowing values 0–9999.  With up to 511
sub-algorithms per `.alg` file and typical move counts of 0–19 per type per
algorithm, no field overflows into the next.

The previous base-1000 encoding (three digits per field) risked overflow when
aggregating costs across many sub-algorithms: at 511 algorithms × up to 20 moves
of one type = 10,220, which exceeds 999.

`avcost` and `wccost` decode with the matching divisors:
```cpp
int v = cost / 1000000000000LL ;
int h = cost / 100000000 % 10000 ;
int o = cost / 10000 % 10000 ;
int c = cost % 10000 ;
```

### Cost-based planning in plan_greedy

`plan_greedy` previously compared candidate plans by **popcount** (number of
algorithms), keeping the plan with fewest activations for each state.  It now
uses the same cost-based comparison as `plan_optimal`:

```cpp
ll nc = curcost[j] + algs[i].cost ;
if (nc < nextcost[d]) {
    nextcost[d] = nc ;
    next[d] = newplan ;
}
```

`curcost[]` and `nextcost[]` vectors track accumulated costs of committed and
candidate plans respectively, reset (`nextcost` and `next`) at the start of
each round.  The commit step remains greedy — only never-set states are
committed:

```cpp
if (!next[j].is_none() && plans[j].is_none()) {
    plans[j] = next[j] ;
    curcost[j] = nextcost[j] ;
}
```

Since `algs[i].cost` encodes `v` at the 10¹² scale (and all current algorithms
have v=2), the primary sort criterion is still algorithm count, with `h` and `o`
as tiebreakers.  The plan stored for each state is therefore the minimum-count
plan with cheapest (h, o) among plans of that count — a more faithful objective
than raw popcount.

`plan_optimal` was already cost-based and required no changes.  Existing
`.prob-greedy` cache files should be deleted when the new planner is first run,
as the plans may differ.

### Memory optimisations: lazy unpack, eager free of plans and unpack

With 28,436 `.alg` files loaded (56,872 netgen objects — one greedy, one
optimal per file), the dominant memory consumers were:

| Structure | Per netgen | Total (all netgens) |
|-----------|-----------|---------------------|
| `plans[]` | goal × 8 B (up to 21.6 MB) | ~133 GB |
| `unpack[]` | goal × 4 B (up to 10.8 MB) | ~67 GB |
| `packhi[]`+`packlo[]` | halfmod × 8 B (≤ 32 KB) | ~300 MB |
| `algs[]`, `stringalgs[]` | STAGE × ~450 B | ~2.8 GB |

The largest goals are C(24,12) = 2,704,156 (`000111`, `010101`) and
C(24,8) = 735,471 (`001111`, `011110`), with 4,620 netgens at the maximum size.

**Three changes eliminate the 200 GB problem:**

**1. `unpack[]` built lazily** (was: eager in constructor).
Previously the constructor filled `unpack[]` alongside `packhi[]`/`packlo[]` in
the same double loop.  With 56,872 netgens constructed during `readalgs()`,
all `unpack` arrays were allocated before any planning — causing 30+ GB of
memory before any `.prob2` file was written.

The fix: `unpack.push_back(...)` and the post-construction consistency check
are removed from the constructor.  `packhi[]`/`packlo[]` are still built there
(they are permanent, needed by `getplan()`).  A new `build_unpack()` method
recomputes `unpack` on demand from `halfmod`, `base`, `basepow.size()`, and
`diggoal`:

```cpp
void build_unpack() {
    const int digval[] = {1, 100, 10000} ;
    int hnx = basepow.size() / 2 ;         // = 2*nonx
    vector<int> digsums(halfmod) ;
    digsums[0] = hnx ;
    for (int i=1; i<halfmod; i++)
        digsums[i] = digval[i % base] + digsums[i / base] - 1 ;
    unpack.reserve(goal) ;
    for (int th=0; th<halfmod; th++)
        for (int bh=0; bh<halfmod; bh++)
            if (digsums[th] + digsums[bh] == diggoal)
                unpack.push_back(th*halfmod+bh) ;
}
```

Both `plan_greedy()` and `plan_optimal()` call `build_unpack()` if
`unpack.empty()` before beginning their BFS loop.  Rebuilding costs
O(halfmod²) — at most 16 M iterations for the C(24,12) groups, taking
milliseconds compared to the full planning time.

`diggoal` is stored as a new `int` member (previously a constructor local);
`nonx` is derived from `basepow.size() / 4` when needed.

**2. `plans[]` and `unpack[]` freed after `ensureprobs()`.**
Once `stageprob` is populated (either read from the `.prob2` cache or computed
via `plan()`), neither `plans` nor `unpack` is needed again until the next
`getplan()` call.  `ensureprobs()` now ends with:

```cpp
plans.clear() ;  plans.shrink_to_fit() ;
unpack.clear() ; unpack.shrink_to_fit() ;
```

**3. `getplan()` made non-const with lazy re-plan guard.**
Removing `plans` from memory created a latent bug: if `getplan()` were called
after `ensureprobs()` (in solving mode), `plans` would be empty and the lookup
would be an out-of-bounds access.  This also affected the warm-cache path,
where `plan()` was never called at all, so `plans` was always empty during
solving.  `getplan()` is now non-const and re-plans on demand:

```cpp
planmask getplan(xorbit &orb) {
    if (plans.empty())
        plan() ;
    ...
}
```

**Result:** loading all 28,436 files now consumes ~3 GB (dominated by `algs[]`
and `stringalgs[]`), down from 30+ GB before any `.prob2` file was written.
After planning completes with a warm cache, memory remains at ~3 GB throughout.

## check_coverage.cpp: Coverage Lower-Bound Prover

`check_coverage.cpp` is a standalone tool that answers the question: *given a
specific candidate pool, is it possible to achieve `max_activations ≤ k` for
any sequence of any depth?*  If not, the tool proves a lower bound
`max_activations ≥ k+1`.

### Mathematical background

Two complementary approaches to prove `max_activations ≥ k+1`:

**Counting argument (any candidate pool, any depth ≤ threshold)**

A sequence of depth D using any candidates covers at most
1 + D + D² + … = (D^(k+1) − 1)/(D − 1) states with ≤ k activations
(each of the D^i states reachable at step i is a unique state under the
best case).  For max_activations ≤ 2 and depth D this is (D² + D + 2)/2.
If this is less than the universe size U, then max_activations ≥ 3 is proven
regardless of which candidates are used.  For C(24,12) = 2,704,156 this holds
for any D ≤ 2325; for C(24,4) = 10,626 it holds for D ≤ 144.

**Set-theoretic (specific candidate pool, any depth)**

Define BFS layers: R₀ = {x₀}, R₁ = {π(x₀) : π ∈ C}, R₂ = {π(s) : π ∈ C, s ∈ R₁}, etc.
Any sequence of depth D with max_activations ≤ k can only reach states in
R₀ ∪ R₁ ∪ … ∪ Rₖ, regardless of how many steps the sequence has.  If
R₀ ∪ … ∪ Rₖ ≠ Universe, then `max_activations ≥ k+1` is proven for **every**
possible sequence of **any** depth that uses only candidates from C.

`check_coverage.cpp` implements this set-theoretic proof for a given k.

### Usage

```
./check_coverage -i INITIAL [-m MAX_PERM_COST] [-k MAX_ACT] [-v] < candidates.costperms
```

- `-i INITIAL` (required) — initial state string, e.g. `00001111111111111111`.
- `-m MAX_PERM_COST` — discard candidates with leading cost > MAX_PERM_COST.
  Use this to restrict the proof to cheap-move pools, making computation tractable.
- `-k MAX_ACT` — max activation depth to test (default: 2).  The tool checks
  whether universe ⊆ R₀ ∪ … ∪ R_k.
- `-v` — verbose: print up to 20 sample unreachable states.
- Reads candidates from stdin in `.costperms` format (same as actbeam10 input).

### Algorithm

1. Read candidates; deduplicate by permutation vector (sort + unique on `p` field).
2. Build the initial state from the `-i` string using the same little-endian
   encoding as actbeam10: `s = Σ sym[i] * base^i`.
3. BFS from R₀ = {x₀}: for each layer Rₗ compute Rₗ₊₁ = {apply(π, s) : π ∈ C, s ∈ Rₗ} \ already covered.
4. Report per-layer counts and wall-clock timing.
5. After k layers: if all universe states are covered, conclude max_activations ≤ k
   is achievable in principle; exit 0.  Otherwise prove max_activations ≥ k+1
   (for this candidate pool, for any depth); exit 1.

**`apply_candidate`** uses little-endian encoding matching actbeam10:
- Base-2: `if ((s >> i) & 1) result |= State(1) << p[i]` (bit-permutation).
- General base: extract symbol `sym = s / base^i % base`, accumulate
  `result += sym * base^{p[i]}`.

### Key finding: proof validity depends on the candidate pool

Whether max_activations ≥ 3 can be proven set-theoretically depends on which
candidates are available.  Two instructive test cases (20-element problem,
universe = C(20,8) = 125,970; **wait** — actually 00001111111111111111 has
n=20, base=2, universe = C(20,4)=4,845 in the tested case):

| Candidate pool | Pool size | R₁ coverage | R₂ coverage | Conclusion |
|----------------|-----------|-------------|-------------|------------|
| Full 10.4M (no cost limit) | 10,402,287 | 53.2% (2577/4845) | 100% (4845/4845) | max_activations ≤ 2 achievable |
| Cost-limited (`-m 2004000`) | 548 | partial | 4060/4845 covered | **PROVED max_activations ≥ 3** |

With the full pool, R2 covers the entire universe — so a solution achieving
max_activations ≤ 2 exists in principle.  With only cheap-move candidates (cost
≤ 2,004,000), 785 states remain unreachable even after two activation layers,
proving no sequence using only those candidates can achieve max ≤ 2.

### Performance note

R₂ computation scales as |C| × |R₁| operations.  With 10.4M candidates and
2577 R₁ states, this is ~26 billion operations — taking ~1475 s on a single CPU
core.  Use `-m` to limit candidates to manageable pools; the proof remains valid
for any sequence that draws only from the restricted pool.

### Relation to search results

When actbeam10/metalbeam1 reports `max_activations 3` with an average very
close to 2, it suggests the worst-case state requires 3 activations but is
rare.  `check_coverage -k 2` with the actual candidate file either:
- confirms the search result is optimal (`max_activations 3` truly necessary), or
- reveals that a depth-2-max solution does exist and the beam search simply
  missed it (increasing `-w`/`-b` may find it).

## Remaining Opportunities

- **Eliminate `best_cost` for uniform-cost runs**: when all candidates have
  cost=1, `best_cost[r] = best_count[r]` always.  Dropping it saves another
  21 MB per `MaterializedState` and eliminates the `snap_cost` copy in phase 2.

- **Separate expansion tool**: a companion program to read `-q` output logs
  and recompute full per-step activation counts and contribution statistics.
  When writing this tool, exploit the near-universal prefix pattern: solutions
  at depth D are almost always the depth D-1 solution plus one candidate.
  Cache the full `compute_activation_counts` DP state (~135 MB: `best_mask`
  86 MB, `best_count` 5 MB, `best_cost` 21 MB, `covered` 338 KB, `states`)
  keyed by `seq`.  Before replaying from scratch, check if any cached `seq` is
  a prefix of the current one; if so, only replay the suffix steps.  At depth D
  with a 1-step suffix this is a D× speedup (e.g. 78× at depth 78).  K=1
  cached entry is almost certainly sufficient; it costs ~135 MB but eliminates
  nearly all replay work across a log file.

- **Modularization**: the shared infrastructure (Codec, Ranker, ChunkLUT, etc.)
  is still duplicated across source files.  See `plans.md` for the full plan.

---

## C(24,4) Binary Coverage: max\_activations ≤ 2

### Problem statement

Switch focus from the C(24,12) general-symbol universe (2,704,156 states) to the
**binary C(24,4)** universe: 24-bit strings with exactly 4 ones, giving
**C(24,4) = 10,626 states**.  Initial state: `000000000000000000001111`.
Candidate pool: `bycosth0cx.costperms` (~80M permutations of 24 elements).
Goal: find a sequence of depth D such that every state is reachable from the
initial state in **at most 2 activations** (i.e., at most 2 candidates applied
in sequence).  Equivalently, R₀ ∪ R₁ ∪ R₂ = Universe where:
- R₀ = {init} (1 state)
- R₁ = {Cᵢ(init) : i ∈ seq} (S1 states)
- R₂ = {Cⱼ(s) : j ∈ seq, s ∈ R₁} (S2 states)

### Theoretical lower bound

The maximum coverage achievable with D candidates and max\_activations ≤ 2 is:
1 + |R₁| + |R₂| ≤ 1 + D + D(D+1)/2

The D(D+1)/2 bound (rather than D²) applies because in practice each ordered
pair (Cᵢ applied first, Cⱼ applied second) tends to produce a state already
covered by the reverse-ordered pair (Cⱼ first, Cᵢ second), halving the effective
distinct coverage.  Setting 1 + D + D(D+1)/2 ≥ 10,626 gives **D ≥ 145** as the
practical minimum.

### Feasibility check: fast\_coverage.cpp

Before searching, confirm that max\_activations ≤ 2 is achievable at all with the
80M-candidate pool.

The original `fast_coverage.cpp` (general-symbol, large `rank_of` table) was
completely rewritten for the binary case:

**Key change — eliminate the 64 MB `rank_of` lookup table.**  The original code
stored `rank_of[s]` for all 24-bit states s (16M entries × 4 bytes = 64 MB).
At R₂ with 80M candidates × 4037 frontier states, this caused catastrophic L3
cache thrashing.  Replaced with a **combinatorial (colex) rank** function using a
BC[25][5] table (500 bytes, fits entirely in L1 cache):

```cpp
static uint32_t BC[25][5];   // binomial coefficients
static inline uint32_t comb_rank(State s) {
    uint32_t r = 0, i = 0, x = s;
    while (x) { int p = __builtin_ctz(x); x &= x-1; r += BC[p][++i]; }
    return r;
}
```

**Results** (4 threads, 80M candidates):
- R₁ (1 step): 4,037 new → 4,037 / 10,626 covered (37.99%) in 0.18 s
- R₂ (2 steps): 6,589 new → 10,626 / 10,626 covered (100%) in 450 s

**Conclusion: max\_activations ≤ 2 is achievable with the full pool.**

### Pool extraction: pool\_extract.cpp

80M candidates is too large to hold in RAM during beam search.  `pool_extract.cpp`
produces a smaller, maximally diverse subset:

- For each candidate C in the pool, compute its **S1 image**: s1 = C(init).
- Keep up to K candidates per distinct s1 rank (default K=100).
- This ensures all 4,037 reachable S1 states are represented.

```
./pool_extract -i 000000000000000000001111 -k 100 < bycosth0cx.costperms \
    > /tmp/pool_24_4_k100.costperms
```

Output: **214,260 candidates** covering all **4,037 distinct S1 states**.
Verified with `fast_coverage`: R2 = 10,626/10,626 in 1.4 s.

### Beam search: cobeam\_b2.cpp

A new beam search program specialized for binary C(n,k) max\_activations ≤ 2.

**Key design decisions vs actbeam10/metalbeam1:**

1. **No ChunkLUT**: for K=4 (4 set bits), applying a permutation to a state costs
   4 ctz+and+shift+or iterations — always faster than building a 64 KB LUT.
   Direct application via `apply_perm(s, p)` is used throughout.

2. **Compact MatState**: each beam state holds only what is needed for k=2 coverage:
   - `best_count[r]` — 0=init, 1=S1, 2=S2, 0xFF=uncovered (10,626 bytes)
   - `covered` — bitset of covered ranks (166 uint64\_t words)
   - `s1_list` — list of S1 states (for S2 expansion)
   - `seq` — candidate indices added so far

3. **Phase 1 scoring (`score_ext`)**: reads parent MatState as const, returns
   `ScoreResult {new_cover, new_s1, sum_delta}` with no side effects.
   Uses a generation-counter Scratch for deduplication within a call.

4. **Phase 2 materialization (`extend_mat`)**: updates MatState in-place.

#### Critical bug fix: incomplete S2 pair tracking

The initial implementation of `extend_mat` only computed **upper-triangular** S2
pairs: applying the new candidate Cₙ to each existing S1 state.  This missed:

- **Lower-triangular**: applying all old candidates C₁…Cₙ₋₁ to the new S1 state
  Cₙ(init)
- **Diagonal**: applying Cₙ itself to Cₙ(init)

Impact: at depth 10, coverage was 56 = 1 + 10 + C(10,2) = 1 + 10 + 45 (only
upper-triangular), instead of the correct 1 + 10 + 100 = 111 maximum.

**Fix**: add `seq` field to `MatState`, push the new candidate index *before*
checking for a new S1 state, then after discovering a new S1 apply all `seq`
candidates (including the new one) to it:

```cpp
m.seq.push_back(ci);
State t1 = apply_perm(m.init_state, p);
uint32_t r1 = comb_rank(t1);
if (m.best_count[r1] > 1) {
    // ... mark t1 as S1 ...
    // Lower-triangular + diagonal: apply all seq candidates to new S1
    for (int cj : m.seq) {
        const uint8_t* pj = perm_of(cj);
        State tnew = apply_perm(t1, pj);
        uint32_t rtnew = comb_rank(tnew);
        if (m.best_count[rtnew] != 0xFF) continue;
        m.best_count[rtnew] = 2; set_covered(m, rtnew);
        ++m.cover_count; m.sum_activations += 2;
    }
}
```

The same fix was applied to `score_ext`, which now also computes lower-triangular
and diagonal S2 gain and tracks exact `sum_delta` (accounting for S2→S1 upgrades
that save 1 activation).

#### Performance notes

- Coverage grows at ~75–80 new states per depth step in the middle range,
  slowing to ~10–15 near saturation.
- Each depth step processes W×B candidate pairs; each `score_ext` call is O(|seq|)
  for the lower-triangular loop.
- `choose_candidates` uses `std::unordered_set<int>` (not `vector<bool>(pool)`)
  to avoid O(pool=214K) zero-initialization per beam state per step.

### Pruning: prune\_b2.cpp

After beam search produces a solution, `prune_b2.cpp` removes redundant candidates:
- Backward pass, forward pass, then N random-order restart passes.
- `compute_coverage` correctly checks all D² S2 pairs (not just upper-triangular).

### Results

| Run | W | B | Seed | Beam depth | Pruned | Time |
|-----|---|---|------|-----------|--------|------|
| fixed1 | 500 | 500 | random | 268 | 267 | 35 s |
| wide2 | 1000 | 1000 | 12345 | 262 | 262 | 214 s |

The pruner removed 0–1 candidates, confirming that every candidate in the beam
search output contributes a **unique S1 state** and is individually necessary.
This is consistent with the lower bound: at D=262, all 262 S1 states are distinct
and the S2 grid covers all 10,625 remaining states.

**Best solution**: depth **262**, coverage **10,626/10,626**, max\_activations **2**,
sum\_activations **20,988** (average 1.975).

Solution log file: `cobeam_b2-262-......-000000000000000000001111.log`
(parseable by `collect_solutions`; shortkey `000001`).

### Format conversion: sol2metalbeam.cpp

`sol2metalbeam.cpp` converts a solution index file + pool into a metalbeam1-format
log with cycle notation and full statistics:

```
./sol2metalbeam -i INITIAL -s SOL_FILE < pool.costperms
```

Outputs: `FOUND depth`, `coverage`, `max_activations`, `sum_sq_activations`,
`sum_activations`, `average_activations`, `sum_activation_cost`,
`average_activation_cost`, `initial`, and one line per candidate in the form:

```
N: candidate IDX cost COST (cycle1)(cycle2)...
```

### Source files

| File | Purpose |
|------|---------|
| `fast_coverage.cpp` | Parallel BFS R1/R2 feasibility check; comb\_rank (no large table) |
| `pool_extract.cpp` | Extract diverse subset: up to K candidates per S1 state |
| `cobeam_b2.cpp` | Beam search for binary C(n,k) max\_activations≤2 |
| `prune_b2.cpp` | Greedy backward/forward/random-restart pruning |
| `sol2metalbeam.cpp` | Convert solution index file to metalbeam1 log format |

---

## metalbeam2.cpp: Improved Scoring and Sort Order

`metalbeam2.cpp` is a fork of `metalbeam1.cpp` incorporating two improvements
derived from experience with `cobeam_b2`:

### Motivation

`metalbeam1` sorts beam states by `sum_sq_activations` (sum of squares) as its
primary criterion after coverage, using `sum_activations` (total/average) as a
secondary tiebreaker.  The sum-of-squares heuristic penalises outliers but does
not directly optimise for the average-case or worst-case experience:

- **Average activations** (`sum_activations / universe`) is the most important
  metric in practice — it determines how often an element is used on a typical
  input.
- **Max activations** is important for worst-case guarantees but is secondary.

Additionally, `score_extension` in `metalbeam1` only scores the **upper-triangular**
S2 contribution of a new candidate: applying the candidate to all existing depth-1
states.  It misses the **lower-triangular** and **diagonal** S2 pairs that arise when
the new candidate creates a new depth-1 state — exactly the bug that was identified
and fixed in `cobeam_b2`.

### Changes

#### 1. Beam sort order (`beam_better`)

Sort order changed from:
```
cover_count↓ → sum_sq_activations↑ → sum_activations↑ → sum_activation_cost↑
```
to:
```
cover_count↓ → sum_activations↑ → max_activations↑ → sum_activation_cost↑
```

`sum_activations` is now the primary optimisation criterion (the average is most
important overall); `max_activations` is secondary (worst case matters, but less
so).  `sum_sq_activations` is no longer used in the sort; it is still computed and
reported for reference.

#### 2. Lower-triangular S2 scoring in `score_extension`

When applying candidate `cand_idx` to the parent `MaterializedState`, if the
application to the init state (the unique depth-0 state) yields a **new** depth-1
state `t1`, `metalbeam2` also applies every candidate in `parent_seq` to `t1` to
discover depth-2 states not yet covered.  The new candidate itself is also applied
to `t1` (the diagonal term).

Implementation:
```cpp
bool     have_new_d1  = false;
State    new_d1_state = 0;
uint64_t new_d1_cost  = 0;

// In the main loop, when new_count == 1 and newly covered:
if (new_count == 1) {
    have_new_d1  = true;
    new_d1_state = s2;
    new_d1_cost  = new_cost;
}

// After the main loop:
if (have_new_d1) {
    auto score_d2 = [&](const ChunkLUT& lut_d2, uint64_t d2_cand_cost) {
        State t2 = apply_lut(new_d1_state, lut_d2, codec);
        uint32_t r2 = ranker.rank(t2);
        uint64_t d2_cost = new_d1_cost + d2_cand_cost;
        // ... score r2 as depth-2 using same temp scratch logic ...
    };
    for (int cj : parent_seq)
        score_d2(lut_cache.get(cj).lut, perms[cj].cost);   // lower-tri
    score_d2(ce.lut, cand_cost);                             // diagonal
}
```

Since permutations are bijective, there is at most one new depth-1 state per
`score_extension` call, so the inner loop over `parent_seq` adds at most
`depth` additional LUT applications — cheap at any practical depth.

This brings `score_extension`'s estimate into agreement with the true score that
Phase 2's `apply_candidate_in_place` would compute, reducing the gap between
beam-selection estimates and final materialised scores.

#### 3. Stopping criterion

The early-exit test is updated to match the new primary sort key:
```cpp
// Before:
if (have_solution && best.sum_sq_activations >= prev_sum_sq_activations)
// After:
if (have_solution && best.sum_activations >= prev_sum_activations)
```

The program stops when average activations stops improving across consecutive
complete solutions, rather than when sum-of-squares stops improving.

### Depth limit raised to 511

The activation mask used in `compute_activation_counts` (the background print
replay) was widened from 256 bits to 512 bits:

```cpp
using Mask256 = std::array<uint64_t, 8>;   // was array<uint64_t, 4>
```

`mask256_zero`, `mask256_popcount`, and `mask256_or` were updated to operate
over all 8 words.  The `-e` limit check and usage string were updated from 256
to 511.

The hot DP path (`score_extension`, `apply_candidate_in_place`) and the Metal
shader never use this mask — they operate only on `best_count`/`best_cost` —
so search speed is unaffected.  The only cost is doubled memory for the two
temporary `best_mask`/`snap_mask` vectors allocated during printing (~172 MB
for n=24, up from ~86 MB), which is acceptable since printing runs in a
background thread.

This allows metalbeam2 to search for solutions of depth up to 511, covering
the C(24,4) solutions found around depth 262.

### Build

```
clang++ -O3 -g -std=c++17 -o metalbeam2 metalbeam2.cpp metal_score.mm \
    -framework Metal -framework Foundation
```

---

## metalbeam3.cpp: Count-Only Beam Search with GPU Scoring

`metalbeam3.cpp` is a new program that searches for the **minimum-depth** sequence
of permutation candidates that achieves full coverage.  It does not track activation
counts or costs at all — the objective is purely to minimise the number of elements
(depth) required.

### Motivation

`metalbeam2` is designed to minimise average/max activations once a full-coverage
sequence exists.  For problems where coverage itself is the only goal (e.g. finding
a short universal traversal before refining activation counts), a count-only searcher
is faster and uses much less memory:

| Per-state field | metalbeam2 | metalbeam3 |
|-----------------|-----------|-----------|
| `best_count`    | 2 B       | **removed** |
| `best_cost`     | 8 B       | **removed** |
| `covered`       | 1/8 B     | 1/8 B |
| `states`        | up to 8 B | up to 8 B |
| **Mat total (T=2.7M)** | **~27 MB** | **~22 MB** |

The `best_count`/`best_cost` arrays (and the `snap_count`/`snap_cost` snapshot
copies in phase 2) are eliminated entirely.  `apply_candidate_in_place` becomes a
simple first-coverage-only loop with no snapshot needed and no improvement branch.

Additionally, the search stops immediately at the first depth that achieves full
coverage rather than continuing to refine activation scores.

### Files

| File | Purpose |
|------|---------|
| `metal_score3.h` | Pure C++ interface; no ObjC headers |
| `metal_score3.mm` | ObjC++ implementation with simplified Metal shader |
| `metalbeam3.cpp` | Count-only beam search with Metal phase 1 |

### Simplified Metal shader (`score_coverage` in `metal_score3.mm`)

The shader is stripped down relative to `beam_score.metal`:

- **Removed**: buffers for `best_count` and `best_cost`; no backwards-pass branch;
  no `tg_sq`, `tg_act`, `tg_cost` reductions.
- **Kept**: `state_of_rank`, `covered`, `lut_data`, `cand_params` (lut_base only,
  no `cand_cost`), rank tables (`high_base`, `high_low_table`, `low_rank_flat`),
  `chunk_mods`, `misc`.
- **Output**: `device int* out_cover` — one `int` per candidate (newly covered delta).
  No `ScoreOut` struct; direct indexed write by `tgid`.

Each threadgroup (one per candidate) strides over all `n_states` with 512 threads,
accumulates a `local_cover` count, reduces via threadgroup shared memory, and writes
the sum.  The CPU then computes `out_cover_count[i] = base_cover_count + raw[i]`.

### `metal_score3.h` / `metal_score3.mm` interface

```cpp
MetalScorer3* metal_scorer3_create(uint32_t n_states,
    const uint64_t* state_of_rank,
    const uint32_t* high_base, const uint32_t* high_low_table, uint64_t high_mod,
    const uint32_t* low_rank_flat, uint64_t low_mod, uint32_t n_low_tables,
    uint64_t split_mod, const uint32_t* chunk_mods, int nchunks);

void metal_scorer3_score(MetalScorer3* ms,
    const uint64_t* covered,
    uint32_t base_cover_count,
    const uint64_t* lut_flat, const uint32_t* lut_offsets,
    int n_candidates,
    uint32_t* out_cover_count);

void metal_scorer3_destroy(MetalScorer3* ms);
```

### metalbeam3.cpp design

Based on `actbeam4` (last version before Mask256) with:

- `BeamState`: `seq`, `cover_count`, `tie`, `parent_idx`, `child_cand`.
  No `sum_activations`, `sum_activation_cost`.
- `MaterializedState`: `covered`, `states`, `cover_count`, `tie`.
  No `best_mask`, `best_count`, `best_cost`.
- `beam_mats[]` maintained across depth steps (metalbeam2 style — no replay
  from scratch).
- **Phase 1 (GPU)**: one Metal dispatch per beam entry scoring all candidate
  branches in parallel.  Returns exact `cover_count` for each child.
- **Phase 1 (CPU fallback)**: multi-threaded scratch-array scoring
  (`score_extension_cpu`) — O(|states|) per call with per-thread `temp_seen[]`.
- **Phase 2**: `apply_candidate_in_place` with no snapshot, no improvement branch
  (first-coverage only).
- `beam_better`: `cover_count > tie` only.
- Stops immediately when `best.cover_count == ranker.total`.

### Correctness of GPU count

Because each candidate permutation is a bijection on the multipermutation state
space, no two currently-covered source states map to the same destination.  The
GPU's per-thread `local_cover += 1` accumulation therefore produces an exact
(not approximate) new-coverage count — no double-counting is possible.

### Measured performance

Test: n=12, base=3, T=34,650 states, 26,149 candidates, `-w 200 -b 200`.

- **GPU confirmed**: `Metal: device="Apple M3 Max"  lut_per_cand=2430`
- **Found depth 19** in 20 s total; each depth step 0.6–1.9 s
  (growing with covered-set size)

### Build

```
clang++ -march=native -O3 -g -std=c++20 metalbeam3.cpp metal_score3.mm \
    -framework Metal -framework Foundation -o metalbeam3
```

---

## brob4.cpp: Recent Planning Improvements

### `-D` flag added to usage

`-D dir` was already implemented (sets `algsdir`, the directory brob4 reads
algorithm files from and writes prob caches to) but was not documented in the
`usage()` output.  Now printed alongside the other options.

### Early termination in `plan_greedy()`

Previously `plan_greedy()` always iterated over all STAGE permutations even
when every goal state had been assigned a plan before the loop completed.  A
`covered` counter (initialised to 1 for the start state) now increments each
time a new state is committed, and the loop exits as soon as `covered == goal`:

```cpp
int covered = 1 ;
for (int i = 0; i < STAGE && covered < goal; i++) {
    ...
    for (int j = 0; j < goal; j++)
        if (!next[j].is_none() && plans[j].is_none()) {
            plans[j] = next[j] ;
            curcost[j] = nextcost[j] ;
            covered++ ;
        }
}
```

No change to `plan_optimal()`, which must examine all permutations to find
globally cheaper plans for already-covered states.

### Used-permutation count in the Use line

`finish_plan()` now counts how many of the STAGE permutations have a non-zero
use count (i.e. appear in at least one state's plan) and stores the result in
a new `int used_perms` field on `netgen` (initialised to -1).

`ensureprobs()` sets `used_perms` from `stageprob` when loading from the prob
cache (where `cnts[]` is not available), using `stageprob[i] > 0` as the proxy.

The `explain()` Use line now appends `used N/M` showing how many of the M
permutations are actually exercised:

```
Use 000111 foo.alg [greedy] cost 1.23 on 000111 map 0 used 47/63
```

If `used_perms` has not been set by the time `explain()` runs (the worst-case
planning path where `ensureprobs()` is never called), `explain()` calls
`plan->net->plan()` first to trigger computation.

### Warning for zero-use permutations in optimal planning

In `finish_plan()`, if `optimal` is true, each permutation `i` with
`cnts[i] == 0` triggers a warning:

```
Warning: optimal plan permutation 7 has zero use count.
```

This indicates a permutation that is never needed in any state's optimal plan —
potentially a sign that the algorithm file contains redundant entries or that
the optimal planner found paths that bypass certain permutations entirely.

---

## beam20 / beam21 / beam22: CPU Beam Search Optimisations

These three programs are progressive performance improvements over `beam19` (the
previous CPU beam search, derived from `actbeam10`), targeting the C(24,12)
problem with n=24, base=2, universe 2,704,156 states.

### beam20: Bitarray-only MaterializedState, No-Copy Scoring

**Key change — `states` vector and `best_cost` removed from `MaterializedState`.**
Coverage is tracked solely as a `std::vector<uint64_t> covered` bitarray plus
`cover_count` and an incremental XOR hash `tie`.  The `states` vector (which
held up to 21 MB of state expansions) is gone.

**`score_candidate` (Phase 1)** reads the parent `MaterializedState` as `const`.
For each rank in the LUT, it tests whether the rank is already covered using
`test_bit(covered, rank)`.  No copy of `MaterializedState` is made.

**Phase 2 parallelised** over the same `-t N` worker threads as Phase 1.
Survivor materialisation is independent per survivor (each writes to its own
`next_mats` slot), so it is safe to parallelise.

**Incremental XOR hash**: `tie ^= hash64(s2)` for each newly-covered rank,
using a 64-bit Fibonacci hash.  This commutative update avoids a full rescan
at each depth step.

Measured speedup: **1.93× over beam19** at `-w 256 -b 256 -t 16`.

### beam21: Inverse Mode at High Coverage

**Key insight**: when more than half the universe is covered, iterating over
*uncovered* states and applying the inverse LUT is cheaper than iterating over
all states.  `beam21` adds a coverage-fraction check to both `score_candidate`
and `apply_candidate_in_place`:

```cpp
if (mat.cover_count * 2 > ranker.total) {
    // inverse mode: iterate uncovered ranks, apply inverse LUT
    auto missing = extract_uncovered_ranks(mat.covered, ranker.total);
    for (uint32_t r : missing) { ... }
} else {
    // forward mode (as in beam20)
    for (uint32_t r = 0; r < ranker.total; ++r) { ... }
}
```

`extract_uncovered_ranks()` scans the bitarray for zero bits; at >50% coverage
this is faster than the forward scan.

Measured speedup: **3.3× over beam19**.

### beam22: Prefetch + Fixed Suffix Branch Scaling

**Software prefetch** added to the hot loops in `score_candidate` and
`apply_candidate_in_place`:

```cpp
static constexpr size_t PF_DIST = 16;
if (i + PF_DIST < n)
    __builtin_prefetch(&ranker.state_of_rank[ranks[i + PF_DIST]], 0, 0);
```

**Suffix search** (`-x N` flag): backward search from prefix using the inverse
LUT on missing states.  At depth `max_depth - suffix_depth`, the search expands
suffix branches and materialises any that achieve full coverage.

Two bug fixes in suffix branch multiplier:

1. **Coverage-percentage threshold**: the 20× branch multiplier (to explore more
   suffix paths) now fires only when `missing < 1%` of universe
   (`start.missing.size() * 100 < ranker.total`), not unconditionally.
   The original unconditional 20× made `-x 2` ~20× slower than `-x 1`
   because at ~80K missing states the expanded search dominated.

2. **`-x 0` accepted**: previously `-x 0` threw an error ("suffix depth must
   be positive").  Now `suffix_depth = 0` is valid and means "no suffix search".
   `final_completion_mode = (suffix_depth > 0)`.

Measured speedup: **~3.5× over beam19**.

### sweep.py: w/b Split Experiment

60 runs comparing (w=512,b=128), (w=256,b=256), (w=128,b=512) at equal budget
w×b = 65,536, 20 runs each.  All 60 runs found depth 27 with nearly identical
costs (~54.12M).  Conclusion: **no measurable w/b distinction** at this budget
level — the problem is too easy for all three configs.

### sweep2.py: Suffix Depth Experiment

Completed experiment: w=8, b=8, d=26, x=0..4, 20 runs each, MAXCOST=4999999.
Results in `sweep2_results.csv`.  Key findings:
- x=0: ~5% success rate (nearly useless)
- x=1: ~55% success rate, ~20 s per run — sweet spot
- x=2: ~60% success rate, ~21 s per run — sweet spot
- x=3: ~20% success rate — counter-productive overhead
- x=4: ~10% success rate — even worse

### `-x 1` Retry: One Extra Forward + Backward Step

When `-x 1` (suffix depth 1) is given and the suffix search at
`depth == max_depth - 1` fails to find a completion, `beam22` now executes one
additional forward step followed by one more suffix search — behaving as if
`-d` were one larger:

**Old behaviour** (suffix fails → break):
```
forward steps 0..max_depth-2
suffix search at max_depth-1 → depth max_depth attempted → NOT FOUND → exit
```

**New behaviour** (suffix fails, suffix_depth==1 → retry once):
```
forward steps 0..max_depth-2
suffix search at max_depth-1 → fails
forward step at max_depth-1           (beam advances one more depth)
suffix search at max_depth   → depth max_depth+1 attempted → report result
```

Implementation: a `bool did_retry_suffix` flag (initially false) gates the
retry.  On the first failure with `suffix_depth == 1`, `max_depth` is
incremented by 1 and `did_retry_suffix` is set to true, then execution falls
through to the normal forward-step code.  On the second failure (or for
`suffix_depth != 1`), the loop breaks normally.

```cpp
if (suffix_depth != 1 || did_retry_suffix) break;
did_retry_suffix = true;
max_depth++;
// fall through to normal forward step, then loop re-fires suffix search
```

Example (`-d 6 -x 1`): forward steps 1–5, suffix search at depth 5 (→ depth 6
attempt, fails), forward step at depth 6, suffix search at depth 6 (→ depth 7
attempt, found with cost 25040014).

---

## collect_solutions.cpp: Recent Updates

### Cost as Fourth Pareto Dimension

`cost` (total algorithm cost) is now a fourth parameter in the Pareto front,
alongside `depth`, `max_activations`, and `sum_activations`.  Lower cost is
better.

**Cost split — base-10,000 fields:**  Each candidate's cost is split into three
fields of 10,000 digits:

```cpp
cur.sum_cost_f1 += cost / 100000000LL;
cur.sum_cost_f2 += (cost / 10000LL) % 10000LL;
cur.sum_cost_f3 += cost % 10000LL;
```

The base was raised from 1,000 (three digits) to 10,000 (four digits) to
prevent field overflow when aggregating costs across many candidates in long
algorithms.  Total cost for Pareto comparison is reconstructed as:

```cpp
sol.total_cost = sol.sum_cost_f1 * 100000000LL
               + sol.sum_cost_f2 * 10000LL
               + sol.sum_cost_f3;
```

The output filename format (`%08lld` for each field) is unchanged.

The updated `dominates()` function:

```cpp
static bool dominates(const Solution& a, const Solution& b) {
    if (a.depth > b.depth || a.max_activations > b.max_activations ||
        a.sum_activations > b.sum_activations || a.total_cost > b.total_cost)
        return false;
    return a.depth < b.depth || a.max_activations < b.max_activations ||
           a.sum_activations < b.sum_activations || a.total_cost < b.total_cost;
}
```

### Skip-if-Existing File Logic

When writing Pareto-front solutions, `collect_solutions` now checks whether each
output file already exists (by `stat()`-ing the path) before writing:

- If the file exists (identified by the FNV-1a hash in the filename), it is
  skipped and counted as **EXISTING**.
- If it is new, it is written and counted as **NEW**.

Each solution line now begins with `NEW` or `EXISTING` instead of `SOLUTION`.
At the end of the run, stderr reports:

```
Written: N new, M already existed
```

This makes repeated runs over growing log collections idempotent and fast.

---

## brob5.cpp: Multi-Threaded Probability Cache Precomputation

`brob5.cpp` is derived from `brob4.cpp` with one structural addition: after
`readalgs()`, all probability caches are pre-computed in parallel using a
thread pool, with a progress report after each file completes.

### Changes

**New includes**: `<thread>`, `<mutex>`, `<atomic>`.

**`io_mutex`** (global `std::mutex`): protects all `cout` output from
`plan_greedy()`, `plan_optimal()`, and the progress reporter.

**Thread-safe timing in plan functions**: `plan_greedy()` and `plan_optimal()`
previously called the global `duration()` (which writes a shared `start`
variable) to measure planning time.  In brob5, each function captures a local
`double t0 = walltime()` at entry and computes elapsed time as
`walltime() - t0` — no shared writes.  The `duration()` call that was in
`plan()` (used only to reset the start marker) is removed entirely.

**`preload_probs()`** (new function, called after `readalgs()`):

```cpp
void preload_probs() {
    vector<netgen*> all_ng;
    for (auto& [key, gens] : algset)
        for (auto* ng : gens)
            all_ng.push_back(ng);
    int total = all_ng.size();
    int nthreads = nprob_threads > 0 ? nprob_threads :
                   thread::hardware_concurrency();
    // ...
    atomic<int> next_idx{0}, completed{0};
    auto worker = [&]() {
        for (;;) {
            int idx = next_idx.fetch_add(1);
            if (idx >= total) break;
            all_ng[idx]->ensureprobs();
            int done = ++completed;
            lock_guard<mutex> lk(io_mutex);
            cout << "  probabilities: " << done << "/" << total << " done" << endl;
        }
    };
    // launch nthreads workers, join all
}
```

Each thread pulls a `netgen*` from an atomic work queue, calls `ensureprobs()`
(which either reads from the `.prob2-greedy`/`.prob2-optimal` cache file or
runs `plan_greedy()`/`plan_optimal()`), then prints progress.  Computations
for different netgens are fully independent (all state is per-netgen).

**`-t N` flag**: sets `nprob_threads` (the number of threads to use for
probability precomputation).  Defaults to `std::thread::hardware_concurrency()`.

**Build:**
```
g++ -march=native -Wall -Wsign-compare -g -std=gnu++20 -O3 -o brob5 brob5.cpp
```

---

## beam22: Additional Infrastructure (Not Previously Documented)

### ThreadLocalLUTCache

Each worker thread owns a `ThreadLocalLUTCache` that maps candidate index → prebuilt
`CachedLUTEntry` (both forward and inverse `ChunkLUT`).  The cache is a fixed-size
slot array indexed by `cand_idx % nslots`; on a miss it rebuilds the LUT from the
permutation image and stores it in the slot (evicting whatever was there).

Slot count is computed at startup from a 1 GB total LUT budget divided equally among
threads, then capped at the candidate pool size:

```cpp
size_t nslots = byte_budget / estimate_lut_entry_bytes(codec);
nslots = std::min(nslots, perms.size());
```

The cache eliminates redundant LUT construction when the same candidate appears in
multiple beam states' branches within a single depth step.

### Cost-Bucket Candidate Sampling

Candidates are grouped by cost into `CostBucket` objects (one bucket per distinct cost
value).  `choose_branch_candidates` first picks a bucket uniformly at random among all
buckets, then picks a candidate uniformly within that bucket.  This gives low-cost
buckets — which are typically much smaller than high-cost ones — a fair chance of
selection rather than being swamped by the large number of expensive candidates.

```cpp
struct CostBucket { uint64_t cost; std::vector<int> candidates; };
// uniform bucket pick, then uniform within-bucket pick
```

The number of distinct cost values is reported as `cost_buckets=N` in the startup
banner.

### Streaming Trim (stream_trim)

Phase 1 child generation uses `maybe_add_child` to prune incrementally rather than
accumulating all W×B candidates and sorting once.  When the local child list reaches a
`trim_threshold` (2× the keep limit), `sort_and_trim_to` runs in-place and the
resulting worst-kept item is cached.  Subsequent candidates worse than that item are
dropped immediately without insertion.  This caps peak memory to ≈2× beam width
regardless of branching factor and avoids the O(W×B log W×B) final sort.

### Pre-extracted Rank Lists

Before Phase 1 scoring, covered and uncovered rank lists are extracted once from each
parent `MaterializedState`:

```cpp
// Phase 1 (scoring)
if (mat.cover_count * 2 > ranker.total)
    uncovered_ranks = extract_uncovered_ranks(mat.covered, ranker.total);
else
    covered_ranks   = extract_covered_ranks(mat.covered, ranker.total);
```

Both `score_candidate` and (after survivor selection) `apply_candidate_in_place`
receive the same pre-extracted lists, avoiding repeated full bitarray scans per
candidate.

### Cheapest-Suffix Mode (cheapest_suffix)

When `-x N` is active, suffix candidates are always drawn from the cheapest available
cost bucket rather than sampled with cost-bucket weighting.  The intent is to keep
suffix sequences short and cheap, since suffix depth is fixed and every step must
contribute coverage.  This is always on (reported as `cheapest_suffix=1` in the banner).

---

## beam24.cpp: Bounded-Activation Beam Search

`beam24.cpp` is derived from `beam22.cpp` and adds a `-k N` option that specifies the
maximum number of activations any state may require.  Under this constraint the search
depth may increase substantially, since coverage is harder to achieve when each state
can only be reached through a chain of at most N candidates applied in sequence.

### Motivation

`beam22` optimises by coverage then cost, with no limit on how many times any
individual candidate must be applied to reach a given state.  In practice many states
require long activation chains.  `beam24` with `-k N` restricts to solutions where
every state is reachable within N applications, enabling quality guarantees at the cost
of deeper searches.

### New global: `max_act_limit`

```cpp
static int max_act_limit = 0;   // 0 = disabled (beam22 behaviour)
```

When 0, beam24 behaves identically to beam22.  When N > 0, the new machinery is
enabled.

### per-state activation count: `best_count`

`MaterializedState` gains a `std::vector<uint8_t> best_count`:

```
0xFF = uncovered
0    = initial state (cost 0 activations)
1..N = covered with that many activations
```

Only allocated when `max_act_limit > 0`.  At full coverage (2.7M states, n=24) this
adds 2.7 MB per beam slot.

`MaterializedState` and `BeamState` also gain `uint64_t sum_sq` — the sum of
`best_count[r]²` over all covered states (excluding the initial state at count=0).

### Secondary optimisation: sum_sq

Beam ordering:
```
cover_count ↓ → sum_sq ↑ → cost ↑ → seq.size() ↑ → tie ↑
```

Minimising the sum of squares naturally penalises states with high activation counts,
driving the worst-case count down without an explicit max term.

### `score_candidate` in max_act mode

Iterates `covered_ranks` (always, in max_act mode — improvement of already-covered
states must be detected).  For each covered rank r0 with `best_count[r0] < N`, applies
the LUT to get r2; if r2 is uncovered or can be improved:

```cpp
uint8_t new_c = c0 + 1;
if (old_c == 0xFF) { ++new_count; delta_sq += new_c * new_c; }
else if (new_c < old_c) { delta_sq += new_c*new_c - old_c*old_c; }
```

Returns a `ScoreResult { cover_count, sum_sq, tie }`.

### Cascade-safe `apply_candidate_in_place`

A subtle correctness issue arises in max_act mode: if candidate C maps r0 → r2 → r3
in two steps, processing r0 first sets `best_count[r2] = c0+1`, then processing r2
again (from the same covered_ranks list) would set `best_count[r3] = c0+2`, effectively
applying C twice from a single call.

Fixed by a two-pass snapshot approach:

1. **Pass 1** — collect `(r0, c0)` pairs from covered_ranks using the current
   (pre-update) `best_count` values.
2. **Pass 2** — apply updates using the snapshotted `c0` values.

This bug does not exist in beam22 because beam22 has no `best_count`; its
`parent_ranks` list is fixed before the loop and contains only pre-existing covered
states.

### Unlimited default depth

The default `max_depth` was changed from 200 to `std::numeric_limits<int>::max()`,
so omitting `-d` truly means unlimited depth.  This is essential for `-k N` runs
which can require hundreds of steps.

### Options

```
-k N    set max activations limit (0 = disabled, default)
```

Banner reports `max_act=N`; progress lines include `sum_sq=N` when max_act is enabled.

### Build

```
g++ -O3 -g -o beam24 beam24.cpp
```

### Results

| k | Width | Branches | Depth found | Notes |
|---|-------|----------|-------------|-------|
| 5 | 40 | 40 | plateau at 2,704,011/2,704,156 | 145 states require >5 activations; wider beam does not help |
| 5 | 80 | 80 | same plateau | confirmed: not a beam-width limitation |
| 6 | 40 | 40 | **FOUND depth 70** | all 2,704,156 states covered, sum_sq=66,288,627 |

The k=5 plateau is a fundamental property of the candidate pool with that activation
limit: the 145 unreachable states require more than 5 activation steps regardless of
sequence length or beam width.

---

## brob15.cpp: Lazy Streaming for Parallel Readmoves

`brob15.cpp` is derived from `brob14.cpp` and adds parallel strip solving when
an input move sequence is supplied via `-i`.  Previously `-i` forced
`alwaysinmem=1` (all center rows allocated up front), making multi-threaded
solving impossible.  Now when `readmoves && fast && nprob_threads > 1`, each
strip is instantiated lazily from the moveset — rows allocated and freed exactly
once, with no allocate-then-free-then-reallocate waste.

### Motivation

For large cubes (N=100, ~49 strips of ~49 rows each), allocating the full 49×49
center row table before solving is prohibitive.  The serial path (1 thread) still
uses `alwaysinmem=1` for simplicity.  The parallel path (>1 thread, fast mode)
uses `alwaysinmem=0` and has each `StripSolver` reconstruct its strip's rows from
the moveset, eliminating peak memory for all strips simultaneously.

### Key changes

**`alwaysinmem` condition** (was: forced on for any readmoves):

```cpp
if (!fast || writemoveinc || randmoves || (readmoves && nprob_threads <= 1))
    alwaysinmem = 1 ;
```

**`g_prep_recorder`** — a global `vector<mv>*` that, when non-null, captures
every move passed to `domove`.  Used during the prep phase of `solve()` to
record the corner-fix moves and wing-cycle moves that modify center rows.

**`domove` modification**: after normalizing the move (dep canonicalization), if
`g_prep_recorder != nullptr`, the normalized move is appended to the vector.

**Prep recording in `solve()`**: two vectors are recorded:

- `prep_before_norm`: face moves from `solvecorners()` (corner fix).
- `prep_after_norm`: moves from the parity fix and `do_wing_cycles()`.

The recorder is active for both phases, separated by the `normalizefacerots()`
call in between.

**`StripSolver::apply_strip_move`** — static method handling all 18 move cases.
Dispatches to `s_rot4hor`/`s_rot2hor`/`s_rot4ver`/`s_rot2ver` (strip-local
variants of the global rotation helpers that take a local `int fr[6]` instead of
`cube.facerots`).  Face moves accumulate into `fr`; slice moves reference `fr`
for position indexing.

**`apply_fr_dosolve` vs `apply_fr_normfr`** — two separate fr-apply helpers with
different rotation conventions, matching the serial path exactly:

| Helper | Convention | Matches |
|--------|------------|---------|
| `apply_fr_dosolve` | `std::rotate` (LEFT rotate by fr) | `dosolve()` |
| `apply_fr_normfr` | `rot4`/`rot2`/`rot4i` (RIGHT rotate by fr) | `normalizefacerots()` |

**Critical bug**: using a single `std::rotate` helper for both the original
moveset and the corner-fix phase gave wrong results for scrambles involving face
moves.  `normalizefacerots` uses `rot4` (right rotation by fr) while `dosolve`
uses `std::rotate` (left rotation by fr) — they are INVERSE operations for the
same fr value.  Using the wrong convention for the corner-fix phase rotated center
stickers in the wrong direction, producing a different (harder-to-solve) center
state that required more moves.

**`fill_strip_from_moveset`** — uses three separate `fr` arrays for the three
phases:

1. `fr_orig`: replay `orig_mvs` (original moveset), then apply `apply_fr_dosolve`
2. `fr_corner`: replay `prep_before` (corner-fix face moves), then apply `apply_fr_normfr`
3. `fr_after`: replay `prep_after` (parity fix + wing cycles); net fr=0, so the
   final `apply_fr_dosolve` is always a no-op

**Modified `solve_strip`**: accepts optional `const vector<mv>*` pointers for
`norm_mvs`, `prep_before`, and `prep_after`.  When `!alwaysinmem && norm_mvs`,
calls `fill_strip_from_moveset` instead of `fill_row`.

**Modified `solve()`** parallel worker: pre-normalizes the moveset (canonicalizes
dep so that dep ≤ Mx+My+3−dep) once before spawning threads.  Each thread's
`solve_strip` call receives the three moveset vectors.

### Correctness verification

- **Simple scramble** ("2U"): 1-thread and 4-thread agree exactly.
- **Complex scramble** ("2U 3R 2F 2L 3B 2D U R F L B D 2U' 3R' 2F' 2L' 3B' 2D'",
  N=100): before the two-convention fix, 1-thread gave 16113 moves while 4-thread
  gave 15075 moves (wrong).  After the fix, both give **15973 moves** exactly.
- 8-thread results also match 4-thread (streaming path is deterministic per strip).

### Build

```
g++ -march=native -Wall -Wsign-compare -g -std=gnu++20 -O3 -o brob15 brob15.cpp
```

---

## collect_solutions.cpp: Error Handling Improvements

### Unknown options exit with usage

Previously, any unrecognised flag (e.g. `-D`) was silently treated as a log filename,
causing confusing behaviour (the filename read would fail and the output directory would
stay at the default `algs3`).  Now the argument parser detects flags:

```cpp
} else if (argv[i][0] == '-') {
    fprintf(stderr, "Unknown option: %s\n", argv[i]);
    usage(argv[0]);
    return 1;
}
```

### Unopenable log files exit immediately

Previously a log file that could not be opened printed a warning and continued.  Now:

```cpp
if (!f) {
    fprintf(stderr, "error: cannot open %s\n", path);
    exit(10);
}
```

Exit code 10 distinguishes this from the normal `return 1` error exits, making it
easier to detect in scripts.

---

## findpermnet.cpp: Beam Search over Costalgs Input

`findpermnet.cpp` is the successor to `beam22.cpp` / `beam25.cpp` for finding
master sequences for brobdicube reduction steps.  The key change from earlier
beam programs is the input format and the orbit-key-based state space.

### Input format: costalgs

`findpermnet` reads lines in **costalgs format**, as produced by
`algutil -c -d 1000 -pc`:

```
COST (24-slot_cycles) algorithm_text
```

For example:
```
2002004000 (0 9 18) U 2R U' 3L U 2R' U' 3L'
```

The cost is a packed integer (base-10000 fields encoding outer-turn, vertical-
slice, and horizontal-slice move counts).  When `-d 1000` is used in algutil,
algorithms involving two depth-3 moves receive a cost with a `2 000 000 000`
prefix that `findpermnet` strips automatically in key mode.

### Two input modes (selected by `-i`)

**Permutation mode** (STRING is not a 6-char all-digit/x string):

- Reads raw permutations and costs.
- STRING is the initial multipermutation (e.g. `11110000`, `001122`).

**Key mode** (STRING is exactly 6 chars, each a digit or `x`/`X`):

- Reads costalgs format.
- The key specifies which of the 24 orbit positions are active (non-`x` faces)
  and how many color groups there are.
- `findpermnet` remaps the full 24-slot permutations to a compact internal
  representation covering only the active positions.
- Strips the `2 000 000 000` cost prefix automatically.
- Examples: `010101` (all faces, 2 color groups), `0x1x2x` (U/F/B only, 3
  color groups).

Key mode detection (line ~2021 in findpermnet.cpp):

```cpp
auto is_key_char = [](char c){ return c == 'x' || c == 'X' || std::isdigit((unsigned char)c); };
bool key_mode = (init.size() == 6 &&
                 std::all_of(init.begin(), init.end(), is_key_char));
```

This correctly identifies `010101` (all digits) as key mode even though it
contains no `x` characters.

### Default thread count

```cpp
int num_threads = (int)std::thread::hardware_concurrency();
```

If `-t N` is not specified, `findpermnet` uses all available hardware threads.

### Output format

Each solution block begins with `initial KEY` and then one line per step:

```
initial 010101
1: cost 2004000 U 2R U' 3L U 2R' U' 3L'
2: cost 2006000 ...
...
```

`pareto` reads this format directly.

---

## pareto.cpp: Pareto Front Extraction

`pareto` replaces `collect_solutions` as the tool that post-processes
`findpermnet` output.  Key differences:

- **No .algs lookup**: reads algorithm text directly from the `N: cost COST text`
  step lines in the `findpermnet` log.  No separate algorithm file needed.
- **Pareto front**: computes a 4-dimensional Pareto front over (depth,
  max_activations, sum_activations, total_cost) — lower is better on all
  dimensions.
- **Output**: writes each Pareto-optimal solution as a `.alg` file to the
  directory specified by `-o DIR`.  Filename encodes key, depth, and cost fields
  so brobdicube can load and compare solutions from multiple runs.
- **Idempotent**: checks for existing files by filename (which encodes content
  via cost fields) before writing, so repeated runs over growing log collections
  are fast.

Usage:
```
pareto -o testalgs permnet-phase1.log permnet-phase2.log
```

---

## algutil.cpp: Embedded 6.tws Puzzle Definition

Previously `algutil` required a `6.tws` file in the working directory (or a
path specified on the command line).  This made it fragile when run from
different directories.

### Change: always use embedded 6.tws

The full content of `6.tws` is now embedded as a raw string literal in
`algutil.cpp`:

```cpp
static const char g_builtin_tws[] = R"TWS(
...puzzle definition...
)TWS";
```

A small edit was made to remove rotation moves that are not needed (saving
~600 lines).

The `read_tws` function signature changed from `read_tws(const string &path)`
to `read_tws(std::istream &f)`, taking a stream instead of a filename.  In
`main`:

```cpp
{ std::istringstream ss(g_builtin_tws); read_tws(ss); }
```

All external file-reading code and any command-line path argument were removed.
Unknown arguments now produce an error and exit rather than being silently
treated as a file path.

`#include <cstdint>` was added (required for `uint64_t` etc. on some platforms).

---

## test_pipeline.py: End-to-End Pipeline

`test_pipeline.py` exercises the full offline pipeline from seed algorithms
through to a brobdicube solve.  The script:

1. Calls `make` to ensure all tools are built.
2. Creates a `test/` subdirectory and `cd`s into it.
3. Runs all pipeline stages, writing intermediate files into `test/`.
4. `cd`s back to the root directory on completion.

### Pipeline stages

| Stage | Command | Output |
|-------|---------|--------|
| Conjugation expansion | `algutil -g -U < ../seed.algs \| algutil -xu` | `expanded.algs` |
| Cost annotation (full) | `algutil -c -d 1000 -pc < expanded.algs` | `expanded.costalgs` |
| Face filter (U/F/B) | `algutil -f 0x0x0x < expanded.algs` | `expanded-0x0x0x.algs` |
| Cost annotation (filtered) | `algutil -c -d 1000 -pc < expanded-0x0x0x.algs` | `expanded-0x0x0x.costalgs` |
| Phase 1 beam search | `findpermnet -b 20 -w 20 -i 010101 < expanded.costalgs` | `permnet-phase1.log` |
| Phase 2 beam search | `findpermnet -b 100 -w 100 -i 0x1x2x < expanded-0x0x0x.costalgs` | `permnet-phase2.log` |
| Pareto extraction | `pareto -o testalgs permnet-phase1.log permnet-phase2.log` | `testalgs/*.alg` |
| Cube solve | `brobdicube -D testalgs -p 512 16384` | stdout |

`seed.algs` must already exist in the `beam/` directory (not auto-generated).

The `-pc` flag in the `algutil -c` commands is critical: it adds cycle notation
to the costalgs output so that `findpermnet` in key mode can parse the
permutation from `(cycles)` fields.

---

## Makefile: Current State

Key changes from earlier versions:

- **`-std=c++20`** (was `gnu++20`): no GNU extensions required.
- **`beam22` and `collect_solutions` removed** from the build (files still
  present in the directory but not built by default).
- **fpng support**: `brobdicube` is compiled with `fpng.cpp` for frame writing.
  The `-w` flag (suppress all warnings) is placed immediately before `fpng.cpp`
  on the command line so it only silences fpng's own warnings (tautological
  comparison, unused function) without suppressing warnings from `brobdicube.cpp`.
- **`clean` removes `test/`** directory created by `test_pipeline.py`.
- **macOS `.dSYM` cleanup**: `clean` removes `.dSYM` bundles on Darwin only,
  using `ifeq ($(UNAME), Darwin)` guard so the Makefile still works on Linux.

```makefile
brobdicube: brobdicube.cpp fpng.cpp
    $(CXX) $(CXXFLAGS) -o $@ brobdicube.cpp -w fpng.cpp $(LDFLAGS)
```

---

## brobdicube.cpp: Linux Portability Fixes

Several GCC (Linux) compilation errors were fixed without affecting clang (macOS)
behavior.

### Member variable shadowing type name

GCC rejects a member variable with the same name as its type.  Two such cases:

1. `struct xcube { threecube threecube; ... }` → renamed member to `tc`:
   - All references `cube.threecube` → `cube.tc`
   - Constructor initializer updated accordingly.

2. `struct Stats { blockalgo *blockalgo; ... }` → renamed member to `ba`:
   - All references `.blockalgo->show()` → `.ba->show()`
   - All `byphase[phasebase].blockalgo` → `byphase[phasebase].ba`

### Missing includes

Added `#include <memory>` (for `std::shared_ptr`) and `#include <array>`.

### Operator precedence bug

```cpp
// Before (bug: & has lower precedence than !=)
if (& mvstat.n != solvecount) { ... }

// After
if (&& mvstat.n != solvecount) { ... }
```

### Dangling pointer warning (GCC only)

`g_prep_recorder = &prep_before_norm` stores the address of a local variable in
a global pointer.  This is safe because `g_prep_recorder` is always cleared to
`nullptr` before the function returns, but GCC's `-Wdangling-pointer` fires.

Fixed with a GCC-only pragma guard:

```cpp
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-pointer"
#endif
   if (streaming) g_prep_recorder = &prep_before_norm ;
   solvecorners(cube) ;
   cube.normalizefacerots() ;
   if (streaming) g_prep_recorder = &prep_after_norm ;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
```

The `!defined(__clang__)` guard is necessary because clang does not know
`-Wdangling-pointer` and will itself warn about the unrecognized pragma.

### Default thread count

```cpp
static int nprob_threads = (int)thread::hardware_concurrency();
```

`brobdicube` now defaults to all hardware threads rather than 1.

---

## fpng: Frame Writing Re-enabled

`fpng.cpp` and `fpng.h` were added to the repository to restore frame-writing
capability in `brobdicube`.  The files are third-party; building them with
`-Wall` produces two harmless warnings (tautological comparison, unused
function).  These are suppressed by placing `-w` immediately before `fpng.cpp`
on the compile command line (see Makefile section above).

---

## pareto.cpp: Corrected v/h/o/c Filename Cost Fields

### Problem

The old three-field cost split (`sum_cost_f1 / sum_cost_f2 / sum_cost_f3`
using base-10⁸ and base-10⁴ boundaries) did not match `algutil`'s encoding,
which uses base-10³ boundaries (`d=1000`): `cost = c + 1000*o + 1000000*h +
1000000000*v`.  The mismatch produced garbage numbers in the output `.alg`
filenames.

### Fix

The three fields were replaced with four properly-named fields:

```cpp
long long sum_v = 0;  // sum of depth-3 inner move counts
long long sum_h = 0;  // sum of depth-2 inner move counts
long long sum_o = 0;  // sum of outer (face) move counts
long long sum_c = 0;  // sum of conjugation wrapper move counts
```

Step-line parsing extracts them with the correct divisors:

```cpp
cur.sum_v += cost / 1000000000LL;
cur.sum_h += (cost / 1000000LL) % 1000LL;
cur.sum_o += (cost / 1000LL) % 1000LL;
cur.sum_c += cost % 1000LL;
```

`total_cost` reconstruction:

```cpp
sol.total_cost = sol.sum_v * 1000000000LL + sol.sum_h * 1000000LL
               + sol.sum_o * 1000LL + sol.sum_c;
```

Output filename now uses four 4-digit fields:

```cpp
snprintf(fname, sizeof(fname), "%s/%s-%04d-%04lld-%04lld-%04lld-%04lld-%s.alg",
         out_dir, sol.shortkey.c_str(), sol.depth,
         sol.sum_v, sol.sum_h, sol.sum_o, sol.sum_c, hash_str);
```

The stdout summary likewise prints `sum_cost=%lld/%lld/%lld/%lld` with
`sum_v, sum_h, sum_o, sum_c`.

---

## collect_solutions.cpp: Same v/h/o/c Fix

The identical four changes (struct fields, `total_cost`, step-line parsing,
filename format) were applied to `collect_solutions.cpp`.  The step-line
format there differs (`" candidate M cost COST"`) but the cost decomposition
is the same.

---

## rename_algs.py: One-Time Bulk Rename Script

`rename_algs.py` renames existing `.alg` files (and any sibling `.prob*`
files) that were written with the old misaligned cost fields to the correct
new four-field `v/h/o/c` format.

### Old vs new filename format

| Format | Pattern |
|--------|---------|
| Old (broken) | `{key}-{depth:04d}-{f1:08d}-{f2:08d}-{f3:08d}-{hash}.alg` |
| New (correct) | `{key}-{depth:04d}-{v:04d}-{h:04d}-{o:04d}-{c:04d}-{hash}.alg` |

### Algorithm

1. Scan the target directory in one `os.listdir` pass, building an `index`
   dict mapping basename (before the first `.`) → list of extensions.
2. Collect all basenames that have a `.alg` extension.
3. Process in batches of 5,000 files per `algutil -c -d 1000` invocation:
   - Concatenate all algorithm lines from the batch into a single stdin block.
   - Decode each output cost as `(v, h, o, c)` with the correct divisors.
   - Sum the per-algorithm costs for each file to get `(sv, sh, so, sc)`.
   - Build the new basename and rename all extensions if the name changed.
4. Fall back to per-file `algutil` calls if a batch fails.
5. Print progress every 1,000 files; final summary of renamed/skipped/errors.

### Options

```
python3 rename_algs.py DIR [--algutil PATH] [--batch N] [--dry-run]
```

`--dry-run` prints what would change without performing any renames.

### Efficiency note

Running `algutil` once per file would be too slow for millions of files.
Batching 5,000 files per call reduces subprocess overhead ~5,000×; the
per-file fallback ensures one bad file cannot block an entire batch.

---

## algutil.cpp: -F Flag (Face-Conjugate Filter)

### Purpose

The `-F` flag drops any input sequence whose first and last tokens are both
pure face moves (alphabetic) and are inverses of each other — that is, the
sequence is "face-conjugated" by a move that cancels from outside.

Example: `U F R U'` is dropped because `U` and `U'` are inverse face moves.

### Implementation

```cpp
static int filter_conj_mode = 0;
// in main loop, after reading line_tokens:
if (filter_conj_mode && line_tokens.size() >= 2) {
    const string &front = key_to_name.at(line_tokens.front());
    const string &back  = key_to_name.at(line_tokens.back());
    if (isalpha((uchar)front[0]) && isalpha((uchar)back[0]) &&
        invert_name(front) == back)
        continue;
}
```

### Bug fix: `line_tokens` not populated without other flags

`line_tokens` was only populated when `expand_mode || cost_mode || g_mode`.
With only `-F` active, the token list was always empty and the filter never
fired.  Fixed by adding `|| filter_conj_mode` to all three occurrences of
that condition (using `replace_all=true`).

---

## algutil.cpp: Default Deduplication Bug Fix

### Problem

The final `else` branch in the output-emit code was:

```cpp
} else {
    if (seen.insert(centers_state).second) { emit...; ++kept; } else { ++dupes; }
}
```

This always deduplicated by centers permutation, even when neither `-u` nor
`-U` was given.  As a side effect, `algutil -c` on a batch was returning fewer
lines than the number of input algorithms — it was silently dropping duplicates.

### Fix

The `else` branch was split into a dedup path (gated on `u_mode`) and a
straight-through path:

```cpp
} else if (u_mode) {
    if (seen.insert(centers_state).second) { emit...; ++kept; } else { ++dupes; }
} else {
    emit_prefixed(pfx, pfxlen, line_start, out_len);
    ++kept;
}
```

Deduplication now only occurs when `-u` or `-U` is explicitly given.

### kept= line suppressed

The diagnostic `fprintf(stderr, "kept=%d dupes=%d errors=%d\n", ...)` line was
commented out (left in place for future debugging).  An `(void)dupes;` was
added to suppress the resulting unused-variable warning.

---

## enumkeys.py: Pipeline Progress, b/w Scaling, Pareto, and brobdicube Steps

### Progress messages in findpermnet loop

Before each `findpermnet` run the script prints:

```
{key}  space={ss:,}  ({nact} active, {ng} groups)  b={b} w={w}  candidates={ncands:,}
```

After the run it parses the log and prints either:

```
  FOUND  depth={d}  cost={c:,}  ({elapsed:.1f}s)
```

or:

```
  NOT FOUND  best={cov:,}/{tot:,}  ({elapsed:.1f}s)
```

### b/w scaling inversely with universe size

All cases now target the same total CPU budget `b * w * universe = 25 * max_ss`
(where `max_ss` is the largest universe among all filtered keys).

Baseline: `b = w = max(5, round(5 * sqrt(max_ss / ss)))`.  This gives b=w=5
for the largest universe and larger values for smaller universes.

When the baseline `b` exceeds the number of available candidates, `b` is capped
at `ncands` and `w` is raised to compensate so the total budget is preserved:

```python
budget = 25 * max_ss
b = max(5, round(5 * (max_ss / ss) ** 0.5))
w = b
if ncands and b > ncands:
    b = ncands
    w = max(b, round(budget / (b * ss)))
```

### Step 4: pareto

After all `findpermnet` runs complete, the script calls:

```python
os.system(f"{P} -o {algs_dir} {' '.join(log_files)}")
```

and reports how many `.alg` files were produced in the output directory.

### Step 5: brobdicube -e

After pareto, the script calls:

```python
os.system(f"{B} -e -D {algs_dir} {CUBE_SIZE}")
```

with `CUBE_SIZE = 16384` to precompute all probability caches, timing and
reporting the elapsed time.

---

## runallreductions.py: Renamed and Extended enumkeys.py

`enumkeys.py` was renamed to `runallreductions.py` and extended with new
options and pipeline improvements.

### New options

**`--dir DIR`** — sets the working directory for all pipeline output files
(default: `pipeline`, resolved relative to the script).  Absolute paths are
used as-is; relative paths are joined to the script's directory.

**`--max-cost N`** — drops algorithms from the costalgs files whose collapsed
cost (with `d=1000`) exceeds N.  Uses the new `algutil -m N` flag (see below)
rather than a separate awk pipe, so filtering happens in a single pass during
cost computation.

### Seed file changed

`bigseed.algs` is now read instead of `seed.algs`.

### brobdicube step: planning sweep instead of fixed size

The final brobdicube step changed from:

```python
os.system(f"{B} -e -D {algs_dir} {CUBE_SIZE}")   # single size precompute
```

to:

```python
os.system(f"{B} -P 4 100000000000000000 -D {algs_dir}")   # sweep 4, 8, 16, ...
```

`-P start end` runs `bestcost` + `explain` for sizes `start, start*2, ...`
up to `end`, giving a cost estimate across the full range of practical cube
sizes rather than precomputing caches for one fixed N.

### Cleanup of expanded* files

After all `findpermnet` runs complete (and before the pareto step), the
script removes the large intermediate files that are no longer needed:

- `expanded.algs`
- `expanded-{fsk}.algs` for every face-subset key
- `expanded-{fsk}.costalgs` for every face-subset key

The pareto step only needs the `permnet-*.log` files, so these can be
deleted immediately.  The script prints a count of files removed.

---

## algutil.cpp: -m N Max-Cost Filter

### Purpose

The `-m N` flag drops any algorithm whose collapsed cost exceeds N.  This
replaces a post-hoc `awk '$1 <= N'` pipe in the pipeline, saving a process
fork and doing the work in the same pass as cost computation.

### Implementation

```cpp
static long long max_cost_filter = -1; // -m N: skip algorithms with collapsed cost > N
```

Arg parsing:

```cpp
else if (string(argv[i]) == "-m") {
    if (++i >= argc) { fprintf(stderr, "-m requires a numeric argument\n"); exit(1); }
    max_cost_filter = atoll(argv[i]);
}
```

`max_cost_filter >= 0` is added to all three `line_tokens` population
conditions (so tokens are populated even without `-c`), and the filter check
runs immediately after cost computation:

```cpp
if (cost_mode || max_cost_filter >= 0) compute_cost(line_tokens, cv, ch, co, cc);
if (max_cost_filter >= 0) {
    long long d = d_multiplier > 0 ? d_multiplier : 1000LL;
    if (collapse_cost(cv, ch, co, cc, d) > max_cost_filter) continue;
}
```

If `-d` is not given, the filter defaults to `d=1000` to match the pipeline
convention.  In practice `-m N` is always paired with `-d 1000`.

---

## brobdicube.cpp: 6-Hyphen Filename Support and Parse Error on Unknown Format

### Bug: new 4-field filenames were silently skipped

The `readalgs()` function distinguished filename formats by counting hyphens
in the stem after the 6-character key:

| Hyphens | Format |
|---------|--------|
| 2 | old: `{key}-{packed_cost}-{hash}.alg` |
| 5 | new: `{key}-{depth}-{f1}-{f2}-{f3}-{hash}.alg` |

After the `pareto.cpp` v/h/o/c fix the filename format became
`{key}-{depth}-{v}-{h}-{o}-{c}-{hash}.alg`, which has **6 hyphens**.
The `else if (nhyphens != 2) continue;` branch silently skipped all new
files, leaving `algset` empty.  `bestcost` then immediately hit the "can't
find any good continuations" error.

### Fix

Accept 5 or 6 hyphens as equivalent depth-based formats (depth is always
the second field at `entry->d_name + 7`):

```cpp
if (nhyphens == 5 || nhyphens == 6) {
    // depth is always the second field (chars 7 onward until next hyphen)
    if (atoi(entry->d_name + 7) > planmask::max_stages())
        continue ;
} else if (nhyphens != 2) {
    fprintf(stderr, "! unrecognised .alg filename format (%d hyphens): %s\n",
            nhyphens, entry->d_name) ;
    error("! unrecognised .alg filename format") ;
}
```

### Error on unrecognised format

Previously any unrecognised hyphen count was silently skipped.  Now it
prints the offending filename and calls `error()` (which aborts), so future
format changes are caught immediately rather than causing mysterious
"no algorithms" failures downstream.

---

## pareto.cpp: -A DIR Directory Analysis Mode

### Purpose

In addition to reading FOUND blocks from log files, `pareto` can now analyse
an existing `.alg` directory and report which files are dominated by others —
without removing anything.

### Usage

```
./pareto -A ALGSDIR
```

### Three-dimensional Pareto front

Each `.alg` file is described by three dimensions (all less-is-better):

| Dimension | Source |
|-----------|--------|
| `depth` | second field of filename (`{key}-{depth:04d}-...`) |
| `total_cost` | collapsed v/h/o/c from filename fields |
| `sum_prob` | sum of per-stage probabilities from `.prob3-greedy` or `.prob3-optimal` |

`sum_prob` (sum of stage activation probabilities) is a proxy for average
activations per orbit and is read from the sidecar probability files.  If
neither `.prob3-greedy` nor `.prob3-optimal` exists for a given `.alg` file
the tool errors out immediately — the probability files are required.

### Output format

For each dominated file:

```
DOMINATED {filename}
  depth=D cost=C sum_prob=P
  by {dominator_filename}
     depth=D cost=C sum_prob=P
```

A summary count line is printed at the end.

### Implementation

New functions added to `pareto.cpp`:

- `DirSol` struct — holds `filename`, `key`, `depth`, `sum_v/h/o/c`,
  `total_cost`, `sum_prob`.
- `parse_alg_filename()` — splits stem on `-` into 7 parts
  (key, depth, v, h, o, c, hash); returns false if not exactly 7 parts.
- `read_sum_prob()` — tries `.prob3-greedy` then `.prob3-optimal`; reads
  `COST PROB` lines, sums the prob column; returns -1.0 if neither file exists.
- `dir_dominates(a, b)` — true if `a.depth ≤ b.depth && a.total_cost ≤
  b.total_cost && a.sum_prob ≤ b.sum_prob` with at least one strict inequality.
- `analyze_algs_dir(dir)` — orchestrates the full analysis.
- `-A ALGSDIR` flag added to arg parsing; calls `analyze_algs_dir()` and returns.

---

## brobdicube.cpp: Robust .alg Filename Error Handling

### Problem: silent skipping of valid .alg files

Two separate bugs caused valid `.alg` files to be silently skipped in
`readalgs()`:

**1. Check order**: The first 6-character key validation and dash check came
*before* the `.alg` extension check.  Non-`.alg` sidecar files (`.prob3-greedy`,
`.eval`, etc.) that happened to start with a valid key prefix would pass the key
check and then hit the extension check — correct behaviour by accident.  But
the extension check should be the *first* filter, with all further checks being
fatal for any file that does pass it.

**2. `atoi` overflow**: For old-format files whose second field is a large packed
cost (e.g. `001xxx-26026052008-...`), `atoi()` produces an implementation-defined
result on overflow.  On this platform `atoi("26026052008")` yielded a positive
value > 511, correctly skipping the file, but `atoi("54054108016")` produced a
negative value, incorrectly passing the depth check and loading the file.
This meant some old-format files were silently dropped while others were
silently accepted — non-deterministic behaviour depending on the packed cost.

### Fix

`readalgs()` restructured:

```cpp
// 1. Extension check first — non-.alg files silently skipped
int len = strlen(entry->d_name);
if (len < 4 || strcmp(entry->d_name + len - 4, ".alg") != 0)
    continue;

// 2. Everything below is an .alg file — format errors are now fatal
int validalgfile = 1;
for (int i = 0; validalgfile && i < 6; i++)
    if (strchr("x0123", entry->d_name[i]) == 0)
        validalgfile = 0;
if (!validalgfile || entry->d_name[6] != '-') {
    fprintf(stderr, "! invalid .alg filename (bad key or missing '-'): %s\n",
            entry->d_name);
    error("! invalid .alg filename");
}

// 3. atoll (not atoi) to avoid 32-bit overflow on large packed-cost fields
if (nhyphens == 5 || nhyphens == 6) {
    long long depth = atoll(entry->d_name + 7);
    if (depth > planmask::max_stages()) {
        fprintf(stderr, "! .alg filename depth %lld exceeds max_stages (%d): %s\n",
                depth, planmask::max_stages(), entry->d_name);
        error("! .alg filename depth exceeds max_stages");
    }
}
```

The result: any `.alg` file that cannot be parsed is a hard error rather than
a silent skip.  This immediately surfaces old-format files (packed-cost in
depth field) that need to be renamed before brobdicube can run.

---

## rename_algs.py: Handle Old Packed-Cost Filenames

### Problem

Some `.alg` files in `algs4` were generated by an older pipeline that encoded
the packed aggregate cost `v*10^9 + h*10^6 + o*10^3 + c` in the second filename
field instead of the depth, yielding names like:

```
001xxx-26026052008-0026-0026-0052-0008-77dd7f67b47c73543112488e1e936127.alg
```

The existing `rename_entry()` function extracted `depth_int = int(parts[1])` =
26026052008 and formatted it with `f'{depth_int:04d}'`, which Python renders as
`26026052008` (`:04d` is a minimum width, not a cap).  The recomputed v/h/o/c
matched the filename fields, so `new_base == base` and the file was silently
skipped with no rename.

A separate class of files had `v=0000` despite using algorithms that have
non-zero v costs — these were created by an old pipeline version that did not
track v separately.  The rename script correctly recomputes v from algutil but
the packed-cost detection was needed first to ensure those renames actually ran.

### Fix

`rename_entry()` gains a `nalgs` parameter (the number of algorithms in the
file) and detects the packed-cost case by threshold:

```python
def rename_entry(base, sv, sh, so, sc, nalgs, exts, algs_dir, dry_run, stats):
    ...
    depth_int = int(depth)
    # Old packed-cost format: field 2 holds v*10^9+h*10^6+o*10^3+c, not depth.
    # The real depth is the number of algorithms in the file.
    if depth_int > 9999:
        depth_int = nalgs
    ...
```

All callers (`process_batch` main path and the per-file fallback) pass `count`
(= `len(algs)`) as `nalgs`.

With this fix, the packed-cost files are renamed to their correct 4-field form:

```
001xxx-26026052008-0026-0026-0052-0008-...  →  001xxx-0013-0026-0026-0052-0008-...
```

And files with v=0000 are corrected when algutil produces a non-zero v sum:

```
x0101x-0462-0000-1674-3621-1948-...  →  x0101x-0462-0924-1674-3621-1948-...
```

Dry-run on `algs4` confirmed 26,895 files will be renamed (783 packed-cost
depth fixes + many v=0 corrections), with no errors.

---

## pareto.cpp: -A Renamed to -D; max_activations Added to Both Modes

### -A renamed to -D

The directory-analysis flag was renamed from `-A ALGSDIR` to `-D ALGSDIR` for
consistency with brobdicube's `-D` flag.  All four occurrences (comment, usage
string ×2, arg-parsing condition) updated.

### Directory mode: max_activations as fourth Pareto dimension

The `.prob3-greedy` / `.prob3-optimal` file format ends with a single-integer
line: `max_activations_seen` (the largest plan size observed during planning,
i.e. the maximum number of stages any orbit state needs).  The old
`read_sum_prob()` function stopped reading at that line (because
`fscanf(f, "%lld %lf", ...)` returns 1 instead of 2 on a single-integer line)
and discarded the value entirely.

The function was renamed `read_prob3()` and extended to capture `max_activations`
from the leftover `cost` variable after the loop exits — exactly the same trick
used internally by `ensureprobs()`:

```cpp
static bool read_prob3(const std::string& dir, const std::string& base,
                       double& sum_prob, int& max_activations) {
    ...
    long long cost = 0;
    while (fscanf(f, "%lld %lf", &cost, &prob) == 2)
        sum += prob;
    fclose(f);
    sum_prob = sum;
    max_activations = (int)cost;   // last line: single integer
    return true;
}
```

`DirSol` gained an `int max_activations` field.  `dir_dominates()` was updated
to a four-dimension check (depth, max_activations, total_cost, sum_prob).
The `DOMINATED` output lines now include `max_act=N`.

### Log-file mode: fix for missing max_activations

In log-file mode, older findpermnet logs may not contain a `max_activations`
line.  The `Solution` struct initialises `max_activations = 0` (the default for
missing fields).  Before this fix, `dominates()` treated 0 as genuinely better
than any real max_activations value, so a solution from an old log could falsely
dominate a solution from a newer log purely because its max_activations appeared
to be 0.

Fix: `dominates()` only uses max_activations as a Pareto dimension when
**both** solutions have a non-zero value:

```cpp
static bool dominates(const Solution& a, const Solution& b) {
    bool both_have_max = (a.max_activations > 0 && b.max_activations > 0);
    if (a.depth > b.depth || a.sum_activations > b.sum_activations || a.total_cost > b.total_cost)
        return false;
    if (both_have_max && a.max_activations > b.max_activations)
        return false;
    return a.depth < b.depth || a.sum_activations < b.sum_activations || a.total_cost < b.total_cost
        || (both_have_max && a.max_activations < b.max_activations);
}
```

---

## pareto.cpp: -k Flag for Undominated-Only Output

### Motivation

When pruning a large algorithm directory, the typical workflow is:
1. Run `pareto -D ALGSDIR` to see which files are dominated.
2. Manually delete the dominated files.

With thousands of files the DOMINATED output is noisy; what's needed is a clean
list of keepers.

### Change

A new `-k` (keep-only) flag suppresses DOMINATED output and instead prints one
line per **undominated** file:

```
filename depth=N max_act=N cost=N sum_prob=N.NNNNNN
```

No summary line is printed.  `-k` is only meaningful with `-D`; in log-file mode
it has no effect.

Implementation: `analyze_algs_dir()` signature changed to
`analyze_algs_dir(const char* dir, bool keep_only)`.  The dominated-output block
is guarded by `if (!keep_only)`; an `else if (keep_only)` branch prints the
one-liner for undominated files.  The final summary block is similarly guarded.
`main()` declares `bool keep_only = false;`, parses `-k`, and passes it through.

---

## pareto.cpp: Dual-Mode Analysis, Effective Depth, Epsilon, and -u Flag

### Per-mode Pareto entries (greedy and optimal treated separately)

Previously each `.alg` file produced one `DirSol` entry, using whichever
`.prob3` sidecar was found first.  Now each file produces **up to two** entries
— one for `[greedy]` and one for `[optimal]` — with `DirSol` gaining a `mode`
field.  A file is reported as dominated only when **all** its entries are
dominated.  Cross-mode domination is allowed: an optimal entry from file B can
dominate a greedy entry from file A, or vice versa.

`read_prob3()` was replaced by `read_prob3_one(mode, ...)` which reads a single
named sidecar.  The loading loop tries both modes and creates a separate `DirSol`
for each found.

### Effective depth replaces nominal depth

The `depth` field in `DirSol` previously came from the second field of the
filename (the nominal algorithm count).  It is now the **count of stages with
non-zero probability** in the prob3 file — i.e. the number of algorithms that
are actually exercised by at least one orbit state's plan.

Under greedy planning (which minimises stage activations), only a small subset
of the nominal algorithm count may be used (e.g. 7 of 49), giving the greedy
entry a much lower effective depth.  Under optimal planning (which minimises
total stage cost), all or most stages tend to be used, giving the optimal entry
a higher effective depth but lower sum_prob.  This means the two entries for the
same file are genuinely non-dominated on the Pareto front.

`read_prob3_one` now counts non-zero-probability lines and writes the result
into `sol.depth`, overriding the value parsed from the filename.

### Epsilon comparison for sum_prob

Stage probabilities are rational numbers with denominators up to ~60 million,
so the smallest meaningful difference between two distinct sum_prob values is
~1/60M ≈ 1.7×10⁻⁸.  Floating-point summation of the same probability values
in different stage orderings can produce sums that differ by a few ULPs
(~10⁻¹⁵), creating spurious Pareto non-domination.

`dir_dominates()` now uses `EPS = 1e-9` when comparing `sum_prob`: a value is
considered "worse" only if it exceeds the other by more than EPS, and
"better" only if it is less by more than EPS.  This collapses FP noise while
preserving real differences.

### -u flag: deduplicate equivalent entries

`-u` (used with `-k -D`) suppresses duplicate output lines where all four
metrics are identical (integers exactly; sum_prob within EPS).  The undominated
entries are sorted by `(depth, max_activations, total_cost, sum_prob)` and only
the first entry in each equivalence class is printed.

`analyze_algs_dir()` signature is now
`analyze_algs_dir(const char* dir, bool keep_only, bool unique_only)`.

---

## brobdicube.cpp: MINWRITE=50 and Loaded-Count in Output

### MINWRITE lowered to 50

`netgen::MINWRITE` was 8192.  Any algorithm set with `goal ≤ 8192` (i.e. fewer
than 8192 orbit states) was re-planned from scratch on every brobdicube run,
since no `.prob3` cache was written.  With MINWRITE=50, caches are written for
all algorithm sets with more than 50 orbit states — effectively all real sets —
so the planning cost is paid only once.

### Algorithm count printed on load

`readalgs()` now tracks how many `.alg` files were successfully loaded and
includes the count in the completion message:

```
Loaded 121867 algs in 29.0s
```

(Previously: `Loaded algs in 29.0s` with no count.)

