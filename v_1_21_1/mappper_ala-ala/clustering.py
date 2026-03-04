import json
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from matplotlib.patches import Ellipse, Rectangle

# =============================================================================
# References for conformation region definitions
# =============================================================================
# [1] Ramachandran, G. N., Ramakrishnan, C., & Sasisekharan, V. (1963).
#     Stereochemistry of polypeptide chain configurations.
#     J. Mol. Biol., 7, 95–99.
#     → Original Ramachandran plot; defines allowed φ/ψ space for peptides.
#
# [2] Zimmerman, S. S., Pottle, M. S., Némethy, G., & Scheraga, H. A. (1977).
#     Conformational analysis of the 20 naturally occurring amino acid residues
#     using ECEPP. Macromolecules, 10(1), 1–9.
#     → Introduced standardized region labels (C5, C7eq, C7ax, alpha_R, alpha_L)
#       based on the type of intramolecular hydrogen-bonded ring or helical sense.
#
# [3] Head-Gordon, T., Head-Gordon, M., Frisch, M. J., Brooks III, C. L.,
#     & Pople, J. A. (1991). Theoretical study of blocked glycine and alanine
#     peptide analogs. J. Am. Chem. Soc., 113(16), 5989–5997.
#     → Ab initio φ/ψ potential surfaces for Ac-Ala-NHMe in vacuo;
#       confirmed C7eq as global minimum and C5 as secondary minimum.
#
# [4] Shi, Z., Olson, C. A., Rose, G. D., Baldwin, R. L., & Kallenbach, N. R.
#     (2002). Polyproline II structure in a sequence of seven alanine residues.
#     Proc. Natl. Acad. Sci. USA, 99(14), 9190–9195.
#     → Canonical PII definition: φ = −75°, ψ = +145° (left-handed helix,
#       3 residues/turn); dominant conformation of alanine peptides in water
#       at low temperature.
#
# [5] Sreerama, N., & Woody, R. W. (2004). On the recognition of finer
#     distinctions in Ramachandran map topology: a study of polyproline II
#     conformation in short alanine peptides. Proc. Natl. Acad. Sci. USA,
#     101(46), 16194–16199.
#     → Confirmed PII center at φ = −75°, ψ = +145°; time-averaged conformation
#       at 40 °C: φ = −80°, ψ = +170°. Supports broader ψ upper bound of ~180°.
#
# [6] Bolhuis, P. G., Dellago, C., & Chandler, D. (2000). Reaction coordinates
#     of biomolecular isomerization. Proc. Natl. Acad. Sci. USA, 97(11),
#     5877–5882.
#     → Transition path sampling of Ac-Ala-NHMe; confirms C7eq and C7ax as
#       distinct metastable states in vacuum; αR minimum in solution at ψ ≈ −20°.
#
# =============================================================================
# Region boundary notes
# =============================================================================
# C5    : Extended β-strand; φ ≈ −155°, ψ ≈ +155°. Five-membered H-bond ring.
#         Refs [2,3]. Boundary: φ ∈ [−180, −130], ψ ∈ [130, 180].
#
# PII   : Polyproline II helix (left-handed); canonical center φ = −75°, ψ = +145°.
#         Refs [4,5]. Updated boundary (this work): φ ∈ [−100, −50], ψ ∈ [120, 180].
#         Upper ψ bound extended to 180° to capture thermally averaged structures
#         (φ ≈ −80°, ψ ≈ +170° at 40 °C, ref [5]).
#
# C7eq  : γ-turn (equatorial methyl); φ ≈ −85°, ψ ≈ +80°. Seven-membered H-bond ring.
#         Refs [2,3]. Boundary: φ ∈ [−110, −50], ψ ∈ [50, 110].
#         NOTE: PII is evaluated before C7eq to avoid misclassification in the
#               overlapping φ ∈ [−100, −50] strip.
#
# alpha_R: Right-handed α-helix; φ ≈ −65°, ψ ≈ −40°.
#         Refs [1,2]. Boundary: φ ∈ [−110, −30], ψ ∈ [−80, −10].
#
# alpha_L: Left-handed α-helix; φ ≈ +60°, ψ ≈ +40°.
#         Ref [1]. Boundary: φ ∈ [30, 90], ψ ∈ [20, 100].
#
# C7ax  : γ-turn (axial methyl); φ ≈ +70°, ψ ≈ −65°. Seven-membered H-bond ring.
#         Refs [2,6]. Boundary: φ ∈ [50, 120], ψ ∈ [−130, −40].

# =============================================================================
# Dihedral angle calculation
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

# =============================================================================
# Region-based conformation assignment (literature-grounded, priority-ordered)
# =============================================================================
def assign_conformation(phi, psi):
    """
    Classify (phi, psi) into a named conformation region.

    Boundaries are based on Zimmerman et al. (1977) [ref 2] and Head-Gordon
    et al. (1991) [ref 3] for C5/C7eq/C7ax/alpha_R/alpha_L, and Shi et al.
    (2002) / Sreerama & Woody (2004) [refs 4,5] for PII.

    Evaluation order matters for overlapping regions:
      PII must precede C7eq (both share φ ∈ [−100, −50]).
    """
    # C5: extended β-strand with five-membered H-bond ring
    if -180 <= phi <= -130 and 130 <= psi <= 180:
        return "C5"
    # C5 periodic boundary wrap (φ ≈ ±180° are equivalent)
    elif -180 <= phi <= -130 and -180 <= psi <= -150:
        return "C5"
    # PII: polyproline II (checked BEFORE C7eq due to φ overlap)
    # Center φ = −75°, ψ = +145° [refs 4,5]; ψ upper bound extended to 180°
    # to capture thermally averaged structures at higher temperature [ref 5].
    elif -100 <= phi <= -50 and 120 <= psi <= 180:
        return "PII"
    # C7eq: γ-turn with equatorial Cβ; seven-membered H-bond ring
    elif -110 <= phi <= -50 and 50 <= psi <= 110:
        return "C7eq"
    # alpha_R: right-handed α-helix
    elif -110 <= phi <= -30 and -80 <= psi <= -10:
        return "alpha_R"
    # alpha_L: left-handed α-helix
    elif 30 <= phi <= 90 and 20 <= psi <= 100:
        return "alpha_L"
    # C7ax: γ-turn with axial Cβ; seven-membered H-bond ring
    elif 50 <= phi <= 120 and -130 <= psi <= -40:
        return "C7ax"
    else:
        return "Other"

# =============================================================================
# Periodic angle wrapping (keeps angles in [−180, 180])
# =============================================================================
def wrap_angle(a):
    return ((a + 180) % 360) - 180

# =============================================================================
# Load reaction network and parse coordinates
# =============================================================================
json_path = "reaction_network.json"
try:
    with open(json_path, "r") as f:
        data = json.load(f)
except FileNotFoundError:
    print("Error: reaction_network.json not found.")
    exit()

phi_angles = []
psi_angles = []
node_ids   = []

for node in data.get("nodes", []):
    xyz_file    = node["xyz_file"]
    local_xyz   = os.path.basename(xyz_file)
    target_file = local_xyz if os.path.exists(local_xyz) else xyz_file

    if not os.path.exists(target_file):
        continue

    try:
        with open(target_file, "r") as f:
            lines = f.readlines()

        coords = [
            list(map(float, line.split()[1:4]))
            for line in lines[2:]
            if len(line.split()) >= 4
        ]
        coords = np.array(coords)

        # Atom indices for Ac-Ala-NHMe (0-indexed, verified against EQ000000.xyz)
        # index 1 → acetyl C=O  (C_prev)
        # index 6 → peptide NH  (N)
        # index 2 → Cα
        # index 4 → peptide C=O (C)
        # index 7 → NHMe N      (N_next)
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

        phi_angles.append(wrap_angle(phi))
        psi_angles.append(wrap_angle(psi))
        node_ids.append(node["node_id"])

    except Exception as e:
        print(f"  Warning: skipping node {node.get('node_id', '?')} — {e}")

if not phi_angles:
    print("Error: No valid coordinates could be parsed.")
    exit()

X = np.column_stack((phi_angles, psi_angles))

# =============================================================================
# Step 1: Region-based assignment (primary classification)
# =============================================================================
region_labels = [assign_conformation(p, s) for p, s in zip(phi_angles, psi_angles)]

# =============================================================================
# Step 2: GMM on "Other" points only (supplementary clustering)
# GMM is restricted to unclassified points to avoid contaminating well-defined
# region assignments with data-driven cluster boundaries.
# Optimal n_components selected by BIC with covariance_type='diag' for
# stability with small sample sizes.
# =============================================================================
other_mask = np.array([lbl == "Other" for lbl in region_labels])
other_X    = X[other_mask]

gmm_labels_other = {}
if other_mask.sum() >= 2:
    best_bic = np.inf
    best_gmm = None
    best_n   = 1
    max_n    = min(6, other_mask.sum())

    for n in range(1, max_n + 1):
        try:
            gmm_candidate = GaussianMixture(
                n_components=n, covariance_type="diag", random_state=42
            )
            gmm_candidate.fit(other_X)
            bic = gmm_candidate.bic(other_X)
            if bic < best_bic:
                best_bic = bic
                best_gmm = gmm_candidate
                best_n   = n
        except Exception:
            pass

    print(f"Optimal GMM components for 'Other' points: {best_n} (BIC={best_bic:.1f})")
    raw_cluster   = best_gmm.predict(other_X)
    other_indices = np.where(other_mask)[0]
    for idx, cluster_id in zip(other_indices, raw_cluster):
        gmm_labels_other[idx] = f"Other_{cluster_id}"

# Merge labels
final_labels = []
for i, base_lbl in enumerate(region_labels):
    if base_lbl == "Other" and i in gmm_labels_other:
        final_labels.append(gmm_labels_other[i])
    else:
        final_labels.append(base_lbl)

# =============================================================================
# Visualization
# =============================================================================
REGION_COLORS = {
    "C5"     : "red",
    "PII"    : "blue",
    "C7eq"   : "magenta",
    "alpha_R": "cyan",
    "alpha_L": "green",
    "C7ax"   : "orange",
}

# Rectangles aligned with updated assign_conformation boundaries
REGION_RECTS = {
    "C5"     : (-180, 130,  50,  50),
    "PII"    : (-100, 120,  50,  60),   # ψ: 120–180, φ: −100–−50  [refs 4,5]
    "C7eq"   : (-110,  50,  60,  60),
    "alpha_R": (-110, -80,  80,  70),
    "alpha_L": (  30,  20,  60,  80),
    "C7ax"   : (  50, -130,  70,  90),
}

fig, ax = plt.subplots(figsize=(10, 10))

# Draw reference region rectangles
for name, (x0, y0, w, h) in REGION_RECTS.items():
    color = REGION_COLORS[name]
    rect  = Rectangle(
        (x0, y0), w, h,
        linewidth=1.5, edgecolor=color, facecolor=color,
        linestyle="--", alpha=0.15
    )
    ax.add_patch(rect)
    ax.text(
        x0 + w / 2, y0 + h + 4, name,
        fontsize=9, ha="center", color=color, fontweight="bold"
    )

# Scatter: one series per unique label for legend
unique_labels = sorted(set(final_labels))
for lbl in unique_labels:
    idx   = [i for i, l in enumerate(final_labels) if l == lbl]
    color = REGION_COLORS.get(lbl, "gray")
    ax.scatter(
        [phi_angles[i] for i in idx],
        [psi_angles[i] for i in idx],
        color=color, label=lbl, s=50,
        edgecolor="black", linewidths=0.6, zorder=4
    )

# Node ID annotations
for i, txt in enumerate(node_ids):
    ax.annotate(
        str(txt),
        (phi_angles[i], psi_angles[i]),
        xytext=(4, 4), textcoords="offset points",
        fontsize=7, color="black"
    )

# Classification summary
print("\n=== Classification Summary ===")
from collections import Counter
counts = Counter(final_labels)
for lbl, cnt in sorted(counts.items()):
    nodes_in_lbl = [str(node_ids[i]) for i, l in enumerate(final_labels) if l == lbl]
    print(f"  {lbl:10s}: {cnt:2d} node(s) → {', '.join(nodes_in_lbl)}")

ax.set_title(
    "Ramachandran Plot: Region-Based Classification with GMM for Outliers\n"
    "(PII boundary updated per Shi et al. 2002 / Sreerama & Woody 2004)",
    fontsize=12
)
ax.set_xlabel(r"$\phi$ (degrees)", fontsize=12)
ax.set_ylabel(r"$\psi$ (degrees)", fontsize=12)
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.grid(True, linestyle="--", alpha=0.4)
ax.axhline(0, color="black", linewidth=0.8)
ax.axvline(0, color="black", linewidth=0.8)
ax.legend(loc="lower right", fontsize=9, framealpha=0.8)

plt.tight_layout()
output_path = "ramachandran_classified.png"
plt.savefig(output_path, dpi=300)
print(f"\nPlot saved as {output_path}")

