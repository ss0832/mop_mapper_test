import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import SpectralClustering

# 1. Load JSON data
with open('reaction_network.json', 'r') as f:
    data = json.load(f)

# Initialize a directed graph to keep forward/reverse barriers
G = nx.DiGraph()
nodes_energy = {}
node_list = []

# Add nodes and store their energies
for node in data['nodes']:
    n_id = node['node_id']
    G.add_node(n_id)
    nodes_energy[n_id] = node['energy_hartree']
    node_list.append(n_id)

# Add forward and reverse edges with respective activation barriers
for edge in data['edges']:
    u = edge['node_id_1']
    v = edge['node_id_2']
    G.add_edge(u, v, weight=edge['barrier_fwd_kcal'])
    G.add_edge(v, u, weight=edge['barrier_rev_kcal'])

# 2. Build Affinity Matrix for Clustering
# We convert the barriers to affinities. Smaller barrier = higher affinity.
n_nodes = len(node_list)
affinity_matrix = np.zeros((n_nodes, n_nodes))
node_to_idx = {n: i for i, n in enumerate(node_list)}

RT_scaling = 5.0
for u, v, d in G.edges(data=True):
    i, j = node_to_idx[u], node_to_idx[v]
    # Use the minimum barrier between two nodes as their undirected connection strength
    barrier = d['weight']
    if G.has_edge(v, u):
        barrier = min(barrier, G[v][u]['weight'])
    
    affinity = np.exp(-barrier / RT_scaling)
    affinity_matrix[i, j] = affinity
    affinity_matrix[j, i] = affinity # Ensure symmetry for clustering

# 3. Perform Spectral Clustering (exactly 3 clusters)
sc = SpectralClustering(n_clusters=3, affinity='precomputed', random_state=42)
cluster_labels = sc.fit_predict(affinity_matrix)

# Map nodes to their assigned clusters
cluster_map = {node_list[i]: cluster_labels[i] for i in range(n_nodes)}

# 4. Define Spatial Layout (3 separate circular rings)
pos = {}
# Define centers for the 3 clusters (forming a large triangle)
centers = [
    np.array([0, 15]),            # Top
    np.array([-13, -7.5]),        # Bottom Left
    np.array([13, -7.5])          # Bottom Right
]

for c_id in range(3):
    # Get nodes belonging to the current cluster
    c_nodes = [n for n in G.nodes() if cluster_map[n] == c_id]
    
    # Create a circular layout for these specific nodes
    sub_G = G.subgraph(c_nodes)
    # Scale adjusts the radius of each cluster's ring
    sub_pos = nx.circular_layout(sub_G, scale=5.0) 
    
    # Shift the circular layout to the cluster's designated center
    for n in c_nodes:
        pos[n] = sub_pos[n] + centers[c_id]

# 5. Visualization setup
fig, ax = plt.subplots(figsize=(20, 16))

# Draw nodes (Maintain energy-based coloring)
energy_values = [nodes_energy[node] for node in G.nodes()]
nx.draw_networkx_nodes(
    G, pos, 
    node_size=800, 
    node_color=energy_values, 
    cmap=plt.cm.coolwarm, 
    alpha=0.9,
    edgecolors='black',
    ax=ax
)

# Draw node labels
nx.draw_networkx_labels(G, pos, font_color='black', font_weight='bold', ax=ax)

# Draw edges and their labels (Maintain curved arrows and numerical barriers)
rad = 0.15
for u, v, d in G.edges(data=True):
    weight = d['weight']
    p1 = np.array(pos[u])
    p2 = np.array(pos[v])
    
    # Draw curved arrow
    ax.annotate(
        "", 
        xy=p2, 
        xytext=p1, 
        arrowprops=dict(
            arrowstyle="->", 
            color="gray", 
            shrinkA=15, 
            shrinkB=15, 
            connectionstyle=f"arc3,rad={rad}",
            linewidth=1.0,
            alpha=0.6
        )
    )
    
    # Approximate placement for the edge label
    diff = p2 - p1
    dist = np.linalg.norm(diff)
    if dist != 0:
        norm = np.array([-diff[1], diff[0]]) / dist
        mid_point = p1 + 0.5 * diff + rad * dist * norm * 0.5
        
        # Add weight text (activation barrier)
        ax.text(
            mid_point[0], mid_point[1], 
            f"{weight:.1f}", 
            fontsize=8, 
            color='darkred',
            ha='center', 
            va='center',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
        )

# Add cluster center labels for clarity
for i, center in enumerate(centers):
    ax.text(
        center[0], center[1], 
        f"Cluster {i+1}", 
        fontsize=24, 
        color='gray', 
        alpha=0.3,
        ha='center', 
        va='center',
        fontweight='bold'
    )

plt.title("Reaction Network: 3 Clusters Grouped Spatially (Node Color = Energy)", fontsize=16)
plt.axis('off')
plt.tight_layout()

# Save and show
plt.savefig("clustered_network_spatial.png", dpi=300, bbox_inches='tight')
plt.show()