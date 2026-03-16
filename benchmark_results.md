# 16-350 HW2 — Motion Planning Report

---

## Compilation

```bash
g++ planner.cpp -o planner.out
```

Example run:
```bash
./planner.out map2.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 1 output.txt
```
---

## Approach

### PRM (Probabilistic Roadmap)

I noticed that even increasing `n` was not enough to help the algorithm explore areas that are hard to get to so I decided to use **Gaussian bridge sampling** to bias node placement near obstacle surfaces. This imrpoved roadmap coverage in narrow passages compared to uniform sampling. Node generation works as follows: two configurations q1 (uniform random) and q2 (Gaussian around q1 with σ = π/10 ≈ 18°) are sampled; a node is added only if exactly one of them is in collision. This concentrates samples near C-space obstacle boundaries where navigation is most constrained.

After generating 2000 nodes, the roadmap connects each node to its k = 160 nearest neighbors (by Euclidean distance in angle space) using linear interpolation at π/20 (9°) resolution as the edge validity check. A* with Euclidean heuristic then finds the shortest path between start and goal through the roadmap.

As a little bonus, I tried using bayesian optimization to find optimal values for the parameters n and k. The results were noisy and after a few trails of running bayesian optimization on 500 samples I found n = 2000 and k = 160 was good. I decided to try using bayesian optimization because I learned the technique through an internship and I wanted to try applying it here.

**Data structures:**
- `Node* graph_set[n]` — flat array of raw `double*` pointers, one per configuration. The start and goal nodes occupy fixed indices 0 and 1 so A* can reference them by index
- `double** graph_matrix` — dense n×n adjacency matrix storing squared Euclidean edge weights. I chose it for O(1) edge lookup during A*; the memory cost (2000² × 8 bytes ≈ 32 MB) is acceptable for n = 2000.
- `std::priority_queue<(distance, index)>` (min-heap) — used during roadmap construction to collect and sort candidate neighbors for each node, so only the k closest are checked for collision
- A* open set: another `std::priority_queue<(f_score, index)>` (min-heap) over node indices, with a parallel `std::vector<double> g_score` and `std::vector<int> parent` for cost and path tracking

### RRT-Connect

The implementation follows the standard bidirectional RRT-Connect algorithm. Two trees are grown simultaneously one from the start configuration and one from the goal and swapped each iteration. On each iteration:

1. A random configuration is sampled uniformly in [0, 2π]^DOF.
2. `Extend` steps the nearer tree toward the sample by `step_size = 1.0` radian (in wrapped angular distance).
3. If the extension succeeds, `Connect` repeatedly extends the other tree toward the new node until it either reaches it (REACHED) or hits an obstacle (TRAPPED).
4. On REACHED, the two partial paths are concatenated and densified at π/20 resolution.

All distance computations use **shortest-arc angular differences** (wrapped to [−π, π]) so that the planner correctly handles angles near the 0/2π boundary.


**Data structures:**
- `struct Node` — holds `std::vector<double> angles` (the joint configuration) and an `int parent` index. Using a vector of angles rather than a raw array helped avoid manual memory management and makes copying nodes safe when the tree vector reallocates
- `std::vector<Node> tree_start` and `std::vector<Node> tree_goal` — each tree is a flat dynamic array of nodes. Parent-child relationships are encoded purely through integer parent indices into the same vector, so path tracing is a simple index walk (`curr = tree[curr].parent`) rather than pointer chasing. This also means the tree structure is never invalidated by reallocations since indices remain stable.
- Two raw `std::vector<Node>*` pointers (`current_tree`, `other_tree`) are swapped each iteration via `std::swap` to alternate which tree is extended, avoiding any copying of the tree data itself.

---

## Benchmark Results

**Map:** `map2.txt` (50 × 50)
**DOFs:** 5
**Methodology:** 5 randomly generated valid start/goal pairs × 4 runs each = 20 runs per planner
**Timeout:** 5 seconds per run

### Start/Goal Pairs

| Pair | Start | Goal |
|------|-------|------|
| 1 | [1.374, 3.175, 0.167, 1.249, 4.083] | [1.463, 0.635, 1.747, 3.994, 2.292] |
| 2 | [0.445, 1.495, 4.203, 1.346, 0.831] | [1.343, 2.520, 0.368, 2.381, 6.191] |
| 3 | [1.122, 6.048, 1.668, 0.681, 2.730] | [1.858, 0.800, 2.642, 5.908, 4.256] |
| 4 | [0.670, 2.671, 1.105, 6.019, 3.254] | [0.487, 1.795, 1.707, 2.009, 3.394] |
| 5 | [1.386, 1.379, 2.738, 0.182, 2.112] | [1.330, 2.054, 4.783, 2.382, 4.725] |

### Summary Statistics

| Metric | PRM | RRT-Connect |
|--------|-----|-------------|
| **(d) Success rate** | 20/20 = **100%** | 20/20 = **100%** |
| **(a) Planning time (ms)** | 434.92 ± 7.17 | **0.10 ± 0.05** |
| **(b) Path cost** | **13.66 ± 7.35** | 15.09 ± 2.83 |
| **(c) Vertices generated** | 2000 ± 0 | **12 ± 6.77** |

*Mean ± standard deviation over successful runs.*

### Per-Run Data — PRM

| Pair | Run | Time (ms) | Cost | Vertices |
|------|-----|-----------|------|----------|
| 1 | 1 | 441.3 | 8.745 | 2000 |
| 1 | 2 | 431.1 | 8.745 | 2000 |
| 1 | 3 | 432.6 | 8.745 | 2000 |
| 1 | 4 | 432.3 | 8.745 | 2000 |
| 2 | 1 | 435.1 | 18.414 | 2000 |
| 2 | 2 | 424.2 | 18.414 | 2000 |
| 2 | 3 | 427.3 | 18.414 | 2000 |
| 2 | 4 | 438.1 | 18.414 | 2000 |
| 3 | 1 | 438.7 | 21.366 | 2000 |
| 3 | 2 | 442.1 | 26.725 | 2000 |
| 3 | 3 | 457.5 | 26.725 | 2000 |
| 3 | 4 | 442.1 | 26.725 | 2000 |
| 4 | 1 | 432.3 | 8.135 | 2000 |
| 4 | 2 | 430.1 | 8.135 | 2000 |
| 4 | 3 | 434.9 | 7.438 | 2000 |
| 4 | 4 | 429.4 | 7.438 | 2000 |
| 5 | 1 | 430.4 | 8.076 | 2000 |
| 5 | 2 | 431.9 | 8.076 | 2000 |
| 5 | 3 | 432.9 | 7.903 | 2000 |
| 5 | 4 | 434.0 | 7.903 | 2000 |

### Per-Run Data — RRT-Connect

| Pair | Run | Time (ms) | Cost | Vertices |
|------|-----|-----------|------|----------|
| 1 | 1 | 0.1 | 16.945 | 7 |
| 1 | 2 | 0.1 | 16.945 | 7 |
| 1 | 3 | 0.1 | 16.945 | 7 |
| 1 | 4 | 0.1 | 16.945 | 7 |
| 2 | 1 | 0.1 | 13.934 | 6 |
| 2 | 2 | 0.1 | 13.934 | 6 |
| 2 | 3 | 0.1 | 13.934 | 6 |
| 2 | 4 | 0.1 | 13.934 | 6 |
| 3 | 1 | 0.2 | 19.408 | 24 |
| 3 | 2 | 0.2 | 19.408 | 24 |
| 3 | 3 | 0.2 | 19.408 | 24 |
| 3 | 4 | 0.2 | 19.408 | 24 |
| 4 | 1 | 0.1 | 11.543 | 9 |
| 4 | 2 | 0.1 | 11.543 | 9 |
| 4 | 3 | 0.1 | 11.543 | 9 |
| 4 | 4 | 0.1 | 11.543 | 9 |
| 5 | 1 | 0.1 | 13.634 | 14 |
| 5 | 2 | 0.1 | 13.634 | 14 |
| 5 | 3 | 0.1 | 13.634 | 14 |
| 5 | 4 | 0.1 | 13.634 | 14 |

---