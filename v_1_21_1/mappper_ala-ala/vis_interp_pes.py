import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import Rbf
from scipy.spatial import cKDTree

# Constants
HARTREE_TO_KCAL_MOL = 627.509

# =============================================================================
# Dihedral angle calculation (Sync with reference implementation)
# =============================================================================
def calc_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle (in degrees) defined by four points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))

def wrap_angle(a):
    """Periodic angle wrapping (keeps angles in [-180, 180])"""
    return ((a + 180) % 360) - 180

# =============================================================================
# Data Extraction
# =============================================================================
def extract_phi_psi(xyz_path):
    """
    Extracts phi and psi dihedrals from an XYZ file using specified indices.
    """
    if not os.path.exists(xyz_path):
        return None, None

    try:
        with open(xyz_path, "r") as f:
            lines = f.readlines()

        coords = [
            list(map(float, line.split()[1:4]))
            for line in lines[2:]
            if len(line.split()) >= 4
        ]
        coords = np.array(coords)

        # 0-indexed definitions
        C_prev_idx = 1
        N_idx      = 6
        Ca_idx     = 2
        C_idx      = 4
        N_next_idx = 7

        phi = calc_dihedral(
            coords[C_prev_idx], coords[N_idx], coords[Ca_idx], coords[C_idx]
        )
        psi = calc_dihedral(
            coords[N_idx], coords[Ca_idx], coords[C_idx], coords[N_next_idx]
        )

        return wrap_angle(phi), wrap_angle(psi)

    except Exception as e:
        print(f"Warning: Failed to parse {xyz_path} - {e}")
        return None, None

# =============================================================================
# Plotting routine
# =============================================================================
def plot_interpolated_pes_masked(eq_points, ts_points, distance_threshold=35.0):
    """
    Generates PES plot using RBF interpolation, masks distant regions,
    and plots EQ/TS nodes along with canonical conformer positions.
    """
    points_for_interp = eq_points + ts_points
    if len(points_for_interp) < 3:
        print("Error: Insufficient data points for interpolation.")
        return

    X, Y, Z = zip(*points_for_interp)
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)

    # Periodic expansion for RBF boundaries (3x3 grid)
    X_ext = np.hstack([X, X-360, X+360, X, X, X-360, X-360, X+360, X+360])
    Y_ext = np.hstack([Y, Y, Y, Y-360, Y+360, Y-360, Y+360, Y-360, Y+360])
    Z_ext = np.hstack([Z] * 9)

    # Grid generation
    grid_phi = np.linspace(-180, 180, 200)
    grid_psi = np.linspace(-180, 180, 200)
    grid_phi_mesh, grid_psi_mesh = np.meshgrid(grid_phi, grid_psi)

    # RBF Interpolation
    try:
        rbf = Rbf(X_ext, Y_ext, Z_ext, function='thin_plate', smooth=0.5)
        final_PES = rbf(grid_phi_mesh, grid_psi_mesh)
    except Exception as e:
        print(f"Error during RBF interpolation: {e}")
        return

    # Distance-based Masking
    grid_points = np.column_stack((grid_phi_mesh.ravel(), grid_psi_mesh.ravel()))
    periodic_data_points = np.column_stack((X_ext, Y_ext))
    tree = cKDTree(periodic_data_points)
    
    distances, _ = tree.query(grid_points)
    distances = distances.reshape(grid_phi_mesh.shape)
    final_PES[distances > distance_threshold] = np.nan

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    
    valid_z = final_PES[~np.isnan(final_PES)]
    if len(valid_z) > 0:
        min_z = np.max([0.0, np.nanmin(valid_z)])
        max_level = np.nanpercentile(valid_z, 95) + 5.0
        levels = np.linspace(min_z, max_level, 40)
    else:
        levels = 40

    contour = ax.contourf(grid_phi, grid_psi, np.clip(final_PES, 0.0, None), 
                          levels=levels, cmap='viridis_r', extend='max')
    
    ax.set_facecolor('whitesmoke')
    cbar = fig.colorbar(contour)
    cbar.set_label('Relative Potential Energy (kcal/mol)', rotation=270, labelpad=20)

    # Plot TS Nodes (Crosses)
    if ts_points:
        ts_X, ts_Y, _ = zip(*ts_points)
        ax.scatter(ts_X, ts_Y, c='black', marker='x', s=35, alpha=0.7, label='TS Nodes', zorder=2)
    
    # Plot EQ Nodes (Hollow circles)
    if eq_points:
        eq_X, eq_Y, _ = zip(*eq_points)
        ax.scatter(eq_X, eq_Y, facecolors='none', edgecolors='black', marker='o', s=45, 
                   linewidths=1.5, alpha=0.9, label='EQ Nodes', zorder=3)

    # Plot canonical conformer positions
    REGION_COLORS = {
        "C5"     : "red",
        "PII"    : "blue",
        "C7eq"   : "magenta",
        "alpha_R": "cyan",
        "alpha_L": "green",
        "C7ax"   : "orange",
    }

    REGION_RECTS = {
        "C5"     : (-180, 130,  50,  50),
        "PII"    : (-100, 120,  50,  60),
        "C7eq"   : (-110,  50,  60,  60),
        "alpha_R": (-110, -80,  80,  70),
        "alpha_L": (  30,  20,  60,  80),
        "C7ax"   : (  50, -130,  70,  90),
    }

    for name, (x0, y0, w, h) in REGION_RECTS.items():
        color = REGION_COLORS.get(name, "black")
        rect  = Rectangle(
            (x0, y0), w, h,
            linewidth=1.5, edgecolor=color, facecolor=color,
            linestyle="--", alpha=0.30, zorder=4
        )
        ax.add_patch(rect)
        ax.text(
            x0 + w / 2, y0 + h + 4, name,
            fontsize=10, ha="center", color=color, fontweight="bold", zorder=5,
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1)
        )
        
    ax.set_xlabel(r'$\phi$ (degrees)')
    ax.set_ylabel(r'$\psi$ (degrees)')
    ax.set_title(f'Interpolated PES with Canonical Conformers (Masked, dist < {distance_threshold}°)')
    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-180, 181, 60))
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Adjust legend to avoid duplicate 'Canonical' labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='lower right', framealpha=0.9)

    fig.tight_layout()
    output_filename = 'interpolated_pes.png'
    fig.savefig(output_filename, dpi=300)
    print(f"Saved PES visualization to `{output_filename}`.")

def main():
    react_net_file = 'reaction_network.json'
    
    if not os.path.exists(react_net_file):
        print(f"Error: File not found -> {react_net_file}")
        return

    with open(react_net_file, 'r') as f:
        data = json.load(f)

    nodes = data.get('nodes', [])
    edges = data.get('edges', [])

    if not nodes:
        print("Error: No nodes found in JSON.")
        return

    min_energy_hartree = min([n.get('energy_hartree', float('inf')) for n in nodes])
    
    eq_points = []
    ts_points = []

    print("Extracting EQ node parameters...")
    for node in nodes:
        if 'xyz_file' in node and 'energy_hartree' in node:
            res = extract_phi_psi(node['xyz_file'])
            if res[0] is not None:
                rel_energy = (node['energy_hartree'] - min_energy_hartree) * HARTREE_TO_KCAL_MOL
                eq_points.append((res[0], res[1], rel_energy))

    print(f"Extracting TS edge parameters (Total edges: {len(edges)})...")
    for edge in edges:
        if 'ts_xyz_file' in edge and 'ts_energy_hartree' in edge:
            res = extract_phi_psi(edge['ts_xyz_file'])
            if res[0] is not None:
                rel_energy = (edge['ts_energy_hartree'] - min_energy_hartree) * HARTREE_TO_KCAL_MOL
                ts_points.append((res[0], res[1], rel_energy))

    print(f"Total points collected: EQ={len(eq_points)}, TS={len(ts_points)}")

    plot_interpolated_pes_masked(eq_points, ts_points)

if __name__ == "__main__":
    main()