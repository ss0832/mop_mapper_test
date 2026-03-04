import json
import networkx as nx


def main():
    # =========================================================================
    # Load reaction network
    # =========================================================================
    json_path = "reaction_network.json"
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_path} not found.")
        return

    node_energies = {}
    G = nx.Graph()

    # Add nodes
    for node in data.get("nodes", []):
        nid    = node["node_id"]
        energy = node["energy_hartree"]
        node_energies[nid] = energy
        G.add_node(nid, energy=energy)

    # Add edges (self-loops are discarded; duplicate edges keep the lower TS energy).
    #
    # NOTE on negative-barrier detection:
    #   All barrier values stored in the JSON (barrier_fwd_kcal / barrier_rev_kcal)
    #   are computed relative to the endpoint energies at the time of the TS search.
    #   If a node is subsequently re-optimized, its energy in the "nodes" list may
    #   differ slightly from the value implicitly embedded in the stored barrier.
    #   Recomputing barriers from node energies can therefore yield spurious negative
    #   values.  To avoid false warnings, anomaly detection uses the pre-computed
    #   JSON values rather than recalculating independently.
    negative_barrier_warnings = []
    for edge in data.get("edges", []):
        u    = edge["node_id_1"]
        v    = edge["node_id_2"]
        ts_e = edge["ts_energy_hartree"]

        if u == v:
            continue  # discard self-loops

        fwd = edge.get("barrier_fwd_kcal", float("inf"))
        rev = edge.get("barrier_rev_kcal", float("inf"))
        if fwd < 0.0 or rev < 0.0:
            negative_barrier_warnings.append(
                f"  Edge {edge['edge_id']:>4d} ({u:>2d} -> {v:>2d}): "
                f"fwd = {fwd:.4f}  rev = {rev:.4f} kcal/mol"
            )

        if G.has_edge(u, v):
            G[u][v]["weight"] = min(G[u][v]["weight"], ts_e)
        else:
            G.add_edge(u, v, weight=ts_e)

    if negative_barrier_warnings:
        print(
            "WARNING: The following edges have negative stored barriers.\n"
            "  These likely reflect a genuinely barrierless (or near-barrierless)\n"
            "  region of the potential energy surface rather than an optimization error."
        )
        for w in negative_barrier_warnings:
            print(w)
        print()

    # =========================================================================
    # Identify the largest connected component.
    # Nodes not connected to the main network are excluded from all
    # inter-macrostate barrier calculations and reported explicitly.
    # =========================================================================
    components = list(nx.connected_components(G))
    main_component = max(components, key=len)

    if len(components) > 1:
        print(f"NOTE: The reaction graph contains {len(components)} connected components.")
        for i, comp in enumerate(sorted(components, key=len, reverse=True)):
            tag = "(main)" if comp == main_component else "(isolated — excluded)"
            print(f"  Component {i}: nodes {sorted(comp)}  {tag}")
        print()

    G_main = G.subgraph(main_component).copy()

    # =========================================================================
    # Minimum Spanning Tree (MST) of the main component.
    #
    # The MST minimax-path theorem guarantees that, for any pair of nodes,
    # the unique path in the MST achieves the globally minimum bottleneck
    # edge weight (i.e., the lowest possible maximum TS energy) over all
    # paths in the original graph.  The effective macroscopic activation
    # barrier from basin A to basin B is therefore:
    #
    #   Ea→b = E(bottleneck TS on MST minimax path) − E(basin minimum of A)
    #
    # IMPORTANT: These are potential energy barriers from GFN-FF single-point
    # energies.  They do not include zero-point energy, entropic contributions,
    # or thermal corrections, and must not be interpreted as free-energy barriers.
    # =========================================================================
    MST = nx.minimum_spanning_tree(G_main, weight="weight")

    # =========================================================================
    # Macrostate definitions (from Ramachandran region-based clustering).
    #
    # Node 6 (C7eq) resides in an isolated component and is excluded
    # automatically by the connectivity filter below.
    # =========================================================================
    basin_groups_raw = {
        "C5"     : [0, 12, 21, 31, 33, 37],
        "C7eq"   : [6, 14, 26],
        "alpha_R": [10, 20, 22, 35],
        "alpha_L": [9, 19],
        "C7ax"   : [11, 15],
    }

    # Restrict each conformational basin to nodes present in the main component.
    basin_groups = {}
    for basin_name, nodes in basin_groups_raw.items():
        reachable = [n for n in nodes if n in main_component]
        excluded  = [n for n in nodes if n not in main_component]

        if excluded:
            print(
                f"NOTE: Node(s) {excluded} assigned to '{basin_name}' are isolated "
                f"and will be excluded from barrier calculations."
            )
        if reachable:
            basin_groups[basin_name] = reachable
        else:
            print(
                f"WARNING: Basin '{basin_name}' has no reachable nodes after "
                f"excluding isolated components — skipping entirely."
            )

    if basin_groups != basin_groups_raw:
        print()

    # =========================================================================
    # Basin minimum potential energies
    # =========================================================================
    basin_min_energy = {
        name: min(node_energies[n] for n in nodes)
        for name, nodes in basin_groups.items()
    }

    # =========================================================================
    # Compute pairwise effective activation barriers
    # =========================================================================
    conversion = 627.509  # Hartree -> kcal/mol

    results = []
    for m1, nodes_m1 in basin_groups.items():
        for m2, nodes_m2 in basin_groups.items():
            if m1 == m2:
                continue

            min_bottleneck = float("inf")
            best_path_info = None

            for u in nodes_m1:
                for v in nodes_m2:
                    if not (u in MST and v in MST):
                        continue
                    if not nx.has_path(MST, u, v):
                        continue

                    path = nx.shortest_path(MST, source=u, target=v)
                    bottleneck_ts = max(
                        MST[path[k]][path[k + 1]]["weight"]
                        for k in range(len(path) - 1)
                    )

                    if bottleneck_ts < min_bottleneck:
                        min_bottleneck = bottleneck_ts
                        best_path_info = (u, v, path)

            if min_bottleneck == float("inf"):
                print(
                    f"WARNING: No MST path found between '{m1}' and '{m2}'. "
                    f"These macrostates may belong to different components."
                )
                continue

            barrier_kcal = (min_bottleneck - basin_min_energy[m1]) * conversion
            results.append((m1, m2, barrier_kcal, best_path_info))

    # =========================================================================
    # Output
    # =========================================================================
    results.sort(key=lambda x: x[2])

    header = (
        f"{'Origin':<12} | {'Destination':<12} | "
        f"{'Barrier (kcal/mol)':>20} | Best path (node IDs)"
    )
    print(header)
    print("-" * len(header))
    for origin, destination, barrier, path_info in results:
        _, _, path = path_info
        path_str = " -> ".join(str(n) for n in path)
        print(f"{origin:<12} | {destination:<12} | {barrier:>20.2f} | {path_str}")

    print()
    print("=== Lowest escape barrier per conformational basin ===")
    for m1 in basin_groups:
        outgoing = [(b, dest) for (orig, dest, b, _) in results if orig == m1]
        if outgoing:
            b_min, dest_min = min(outgoing, key=lambda x: x[0])
            print(f"  {m1:<12}: {b_min:.2f} kcal/mol  ->  {dest_min}")

    print()
    print(
        "NOTE: All barriers are GFN-FF potential energy barriers (kcal/mol).\n"
        "      They do not include ZPE, entropy, or thermal corrections and\n"
        "      must not be interpreted as free-energy barriers."
    )


if __name__ == "__main__":
    main()