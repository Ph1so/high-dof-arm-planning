#!/usr/bin/env python3
"""
Benchmark RRT-Connect and PRM planners on map2.txt.

Generates 5 random valid start/goal pairs, runs each 4 times per planner
(20 total runs per planner).  Reports mean ± std for:
  (a) Planning time
  (b) Path cost
  (c) Vertices generated
  (d) Success rate within 5 seconds
"""

import subprocess
import os
import re
import math
import random
import tempfile
import sys
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SCRIPT_DIR     = Path(__file__).parent
PLANNER_BIN    = SCRIPT_DIR / "planner.out"
MAP_FILE       = str(SCRIPT_DIR / "map2.txt")
NUM_DOFS       = 5
LINKLENGTH     = 10          # LINKLENGTH_CELLS in C++
PI             = math.pi
PLANNER_PRM    = 3
PLANNER_RRT    = 1
N_PAIRS        = 5
N_RUNS         = 4
TIMEOUT_SEC    = 5.0

# ---------------------------------------------------------------------------
# Map loading  (replicates C++ indexing: map[y + x*width])
# ---------------------------------------------------------------------------
def load_map(filepath):
    with open(filepath) as f:
        text = f.read()

    height = int(re.search(r'height\s+(\d+)', text).group(1))
    width  = int(re.search(r'width\s+(\d+)',  text).group(1))

    # Strip header to get raw grid characters
    data  = re.sub(r'height\s+\d+', '', text)
    data  = re.sub(r'width\s+\d+',  '', data)
    chars = [c for c in data if not c.isspace()]

    # Match C++ storage: map[y + x*width]
    grid = [0] * (height * width)
    for y in range(height):
        for x in range(width):
            file_idx = y * width + x
            if file_idx < len(chars) and chars[file_idx] != '0':
                grid[y + x * width] = 1

    return grid, width, height


# ---------------------------------------------------------------------------
# Arm collision checker  (replicates IsValidArmConfiguration)
# ---------------------------------------------------------------------------
def _cell(x, y, xs, ys):
    return max(0, min(int(x), xs - 1)), max(0, min(int(y), ys - 1))


def _line_free(x0, y0, x1, y1, grid, xs, ys):
    if not (0 <= x0 < xs and 0 <= x1 < xs and 0 <= y0 < ys and 0 <= y1 < ys):
        return False
    cx0, cy0 = _cell(x0, y0, xs, ys)
    cx1, cy1 = _cell(x1, y1, xs, ys)
    n = max(abs(cx1 - cx0), abs(cy1 - cy0)) + 1
    for i in range(n):
        t  = i / (n - 1) if n > 1 else 0
        nx = int(round(cx0 + t * (cx1 - cx0)))
        ny = int(round(cy0 + t * (cy1 - cy0)))
        if grid[ny * xs + nx] == 1:   # GETMAPINDEX(nx,ny,xs,ys) = ny*xs + nx
            return False
    return True


def valid_config(angles, grid, xs, ys):
    x1, y1 = xs / 2.0, 0.0
    for a in angles:
        x0, y0 = x1, y1
        x1 = x0 + LINKLENGTH * math.cos(2 * PI - a)
        y1 = y0 - LINKLENGTH * math.sin(2 * PI - a)
        if not _line_free(x0, y0, x1, y1, grid, xs, ys):
            return False
    return True


def random_valid(grid, xs, ys, rng):
    while True:
        angles = [rng.uniform(0, 2 * PI) for _ in range(NUM_DOFS)]
        if valid_config(angles, grid, xs, ys):
            return angles


# ---------------------------------------------------------------------------
# Run planner subprocess
# ---------------------------------------------------------------------------
def run_planner(planner_id, start, goal, out_file):
    """Returns dict: success, time_ms, cost, vertices  (values may be None on failure)."""
    cmd = [
        str(PLANNER_BIN),
        MAP_FILE,
        str(NUM_DOFS),
        ",".join(f"{a:.8f}" for a in start),
        ",".join(f"{a:.8f}" for a in goal),
        str(planner_id),
        out_file,
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT_SEC
        )
    except subprocess.TimeoutExpired:
        return dict(success=False, time_ms=None, cost=None, vertices=None,
                    note="timeout")

    out = proc.stdout + proc.stderr

    if proc.returncode != 0 or "no_path_found" in out:
        vert = re.search(r"vertices:\s*(\d+)", out)
        return dict(success=False, time_ms=None, cost=None,
                    vertices=int(vert.group(1)) if vert else None,
                    note="no_path_found")

    t = re.search(r"time:\s*([\d.]+)\s*ms",  out)
    c = re.search(r"cost:\s*([\d.]+)",        out)
    v = re.search(r"vertices:\s*(\d+)",       out)

    return dict(
        success  = True,
        time_ms  = float(t.group(1)) if t else None,
        cost     = float(c.group(1)) if c else None,
        vertices = int(v.group(1))   if v else None,
        note     = "",
    )


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def _stats(values):
    """Return (mean, std) or (None, None) if no data."""
    vals = [v for v in values if v is not None]
    if not vals:
        return None, None
    arr  = np.array(vals, dtype=float)
    mean = float(np.mean(arr))
    std  = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
    return mean, std


def print_report(results, planner_name):
    total     = len(results)
    successes = [r for r in results if r["success"]]
    n_ok      = len(successes)

    print(f"\n{'='*62}")
    print(f"  {planner_name}   ({total} runs, {n_ok} successful)")
    print(f"{'='*62}")

    # (d) success rate
    print(f"  (d) Success rate  : {n_ok}/{total} = {n_ok/total:.1%}")

    # (a–c) stats over successful runs only
    for label, key, unit in [
        ("(a) Planning time ", "time_ms",  " ms"),
        ("(b) Path cost     ", "cost",     "   "),
        ("(c) Vertices      ", "vertices", "   "),
    ]:
        vals        = [r[key] for r in successes]
        mean, std   = _stats(vals)
        n_valid     = sum(1 for v in vals if v is not None)
        if mean is not None:
            print(f"  {label}: {mean:>12.4f}{unit}  ±  {std:.4f}  (n={n_valid})")
        else:
            print(f"  {label}: no data")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # --- Compile if needed ---
    if not PLANNER_BIN.exists():
        print("Compiling planner.cpp …")
        ret = subprocess.run(
            ["g++", "-O2", "-std=c++17", "-o", str(PLANNER_BIN),
             str(SCRIPT_DIR / "planner.cpp"), "-lm"],
            capture_output=True, text=True,
        )
        if ret.returncode != 0:
            print("Compilation failed:\n", ret.stderr)
            sys.exit(1)
        print("Compilation successful.\n")

    # --- Load map ---
    print(f"Loading map: {MAP_FILE}")
    grid, xs, ys = load_map(MAP_FILE)
    print(f"  Map size: {xs} x {ys}\n")

    # --- Generate valid pairs ---
    rng   = random.Random(42)
    pairs = []
    print("Generating 5 valid start/goal pairs …")
    for i in range(N_PAIRS):
        start = random_valid(grid, xs, ys, rng)
        goal  = random_valid(grid, xs, ys, rng)
        pairs.append((start, goal))
        s_str = "[" + ", ".join(f"{a:.3f}" for a in start) + "]"
        g_str = "[" + ", ".join(f"{a:.3f}" for a in goal)  + "]"
        print(f"  Pair {i+1}  start: {s_str}")
        print(f"          goal:  {g_str}")

    # --- Run planners ---
    planners = [(PLANNER_PRM, "PRM"), (PLANNER_RRT, "RRT-Connect")]
    all_results = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        for planner_id, name in planners:
            print(f"\n{'─'*62}")
            print(f"  Running {name} …")
            print(f"{'─'*62}")
            results = []

            for pi, (start, goal) in enumerate(pairs):
                for ri in range(N_RUNS):
                    out_f = os.path.join(tmpdir, f"{name}_{pi}_{ri}.txt")
                    print(f"  Pair {pi+1} run {ri+1}/{N_RUNS} … ", end="", flush=True)
                    r = run_planner(planner_id, start, goal, out_f)
                    if r["success"]:
                        print(f"OK  time={r['time_ms']:.1f}ms  "
                              f"cost={r['cost']:.3f}  verts={r['vertices']}")
                    else:
                        print(f"FAILED ({r.get('note','?')})")
                    results.append(r)

            all_results[name] = results

    # --- Summary ---
    print(f"\n\n{'='*62}")
    print("  BENCHMARK SUMMARY  (map2.txt, 5 pairs × 4 runs each)")
    print(f"{'='*62}")
    for _, name in planners:
        print_report(all_results[name], name)
    print()


if __name__ == "__main__":
    main()
