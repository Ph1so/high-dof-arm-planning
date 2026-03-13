"""
Graph Visualizer for Robotic Arm Planning (PRM and RRT-Connect).

Usage:
    python3 visualize_prm.py <map_file> [graph_file]

    map_file   - the obstacle map (e.g. map1.txt)
    graph_file - the graph/tree file written by the planner.
                 Defaults to prm_graph_<mapname>.txt if not specified.

Examples:
    python3 visualize_prm.py map1.txt                          # PRM
    python3 visualize_prm.py map1.txt rrt_connect_graph_map1.txt  # RRT-Connect

The graph file format (DOFS/NODES/EDGES) is the same for both planners.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc

LINK_LENGTH = 10  # cells, must match LINKLENGTH_CELLS in planner.cpp
PI = np.pi


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_map(map_file: str):
    """Return (grid, width, height).
    grid[y][x] == 1 means obstacle, 0 means free.
    y=0 is the top row (matches image/GIF convention).
    """
    with open(map_file, "r") as f:
        height = int(f.readline().split()[1])
        width  = int(f.readline().split()[1])
        grid = []
        for line in f:
            vals = [int(v) for v in line.split()]
            if vals:
                grid.append(vals)
    return np.array(grid, dtype=np.uint8), width, height


def load_prm_graph(prm_file: str):
    """Return (nodes_dict, edges_list, dofs).
    nodes_dict: {node_index: [angle0, angle1, ...]}
    edges_list: [(i, j), ...]
    """
    nodes = {}
    edges = []
    dofs  = 0
    mode  = None

    with open(prm_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("DOFS"):
                dofs = int(line.split()[1])
            elif line.startswith("NODES"):
                mode = "nodes"
            elif line.startswith("EDGES"):
                mode = "edges"
            elif mode == "nodes":
                parts = line.split()
                nodes[int(parts[0])] = [float(v) for v in parts[1:]]
            elif mode == "edges":
                parts = line.split()
                if len(parts) == 2:
                    edges.append((int(parts[0]), int(parts[1])))

    return nodes, edges, dofs


# ---------------------------------------------------------------------------
# Arm geometry  (mirrors IsValidArmConfiguration in planner.cpp)
# ---------------------------------------------------------------------------

def arm_joints(angles, x_size):
    """Return list of (x, y) joint positions including the base.

    Matches the C++ forward kinematics:
        x1 = x0 + LINKLENGTH * cos(2*PI - theta)
        y1 = y0 - LINKLENGTH * sin(2*PI - theta)
    which simplifies to:
        x1 = x0 + L * cos(theta)     [cos(2π-θ) = cos θ]
        y1 = y0 + L * sin(theta)     [-sin(2π-θ) = sin θ, y increases downward]
    """
    pts = [(x_size / 2.0, 0.0)]
    x0, y0 = pts[0]
    for theta in angles:
        x1 = x0 + LINK_LENGTH * np.cos(theta)
        y1 = y0 + LINK_LENGTH * np.sin(theta)
        pts.append((x1, y1))
        x0, y0 = x1, y1
    return pts


# ---------------------------------------------------------------------------
# Main visualizer
# ---------------------------------------------------------------------------

def visualize(map_file: str, graph_file: str = None):
    basename = os.path.basename(map_file)

    if graph_file is None:
        graph_file = f"prm_graph_{basename}"

    if not os.path.exists(graph_file):
        print(f"ERROR: Graph file not found: {graph_file}")
        print("Run the planner first to generate the graph file.")
        sys.exit(1)

    graph_basename = os.path.basename(graph_file)
    planner_label = graph_basename.replace(f"_{basename}", "").replace("_graph", "").replace("_", " ").strip().upper()
    if not planner_label:
        planner_label = "Planner"

    grid, width, height = load_map(map_file)
    nodes, edges, dofs = load_prm_graph(graph_file)

    if not nodes:
        print("No nodes in graph.")
        sys.exit(1)

    print(f"Map: {width}x{height}  |  Nodes: {len(nodes)}  |  Edges: {len(edges)}  |  DOFs: {dofs}")

    # -----------------------------------------------------------------------
    # Figure: map as background
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 8))

    # The C++ loadMap stores data as map[y + x*width], so file rows = physical x
    # and file columns = physical y. Transpose so imshow sees grid[y][x].
    ax.imshow(grid.T, cmap="binary", origin="upper",
              extent=[0, width, height, 0], zorder=0)

    # -----------------------------------------------------------------------
    # Draw every arm configuration as ghost lines
    # -----------------------------------------------------------------------
    # Collect all link segments across all nodes into one LineCollection
    # for efficient rendering.
    all_segments = []
    for angles in nodes.values():
        pts = arm_joints(angles, width)
        for k in range(len(pts) - 1):
            all_segments.append([pts[k], pts[k + 1]])

    lc = mc.LineCollection(
        all_segments,
        colors=(0.0, 0.6, 0.3),   # dark green, like the GIF
        linewidths=0.6,
        alpha=0.12,                # very transparent — many arms overlap
        zorder=1,
    )
    ax.add_collection(lc)

    # -----------------------------------------------------------------------
    # Highlight start (node 0) and goal (node 1) in bold
    # -----------------------------------------------------------------------
    for idx, color, label in [(0, "lime", "Start"), (1, "red", "Goal")]:
        if idx in nodes:
            pts = arm_joints(nodes[idx], width)
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.plot(xs, ys, color=color, linewidth=2.0, zorder=3, label=label)
            ax.plot(xs[0], ys[0], "o", color=color, markersize=6, zorder=4)

    # -----------------------------------------------------------------------
    # Cosmetics
    # -----------------------------------------------------------------------
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)   # y=0 at top
    ax.set_xlabel("x (cells)", fontsize=10)
    ax.set_ylabel("y (cells)", fontsize=10)
    ax.set_title(
        f"{planner_label} Sampled Configurations — {basename}\n"
        f"{len(nodes)} nodes · {len(edges)} edges · {dofs} DOFs",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="lower right")

    stem    = os.path.splitext(graph_basename)[0]
    out_png = f"visual_{stem}.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved to '{out_png}'")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 visualize_prm.py <map_file> [graph_file]")
        print("  e.g. python3 visualize_prm.py map1.txt")
        print("  e.g. python3 visualize_prm.py map1.txt rrt_connect_graph_map1.txt")
        sys.exit(1)
    graph_arg = sys.argv[2] if len(sys.argv) >= 3 else None
    visualize(sys.argv[1], graph_arg)
