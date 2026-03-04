Although the exploration was configured for 1,000 iterations, the execution was prematurely terminated at iteration 411 following an unexpected operating system failure.

<memo>

To use GFN-FF,
```
pip install pygfnff
```
and
```
touch software_path.conf
```

To execute the reaction pathway exploration, the configuration file (config_snapshot.json) and the input coordinate file (ala_dipeptide_gfnff.xyz) were placed within the same working directory. The tool was invoked using the run_mapper.py script provided by the MultiOptPy package.

Command:
```
python run_mapper.py ala_dipeptide_gfnff.xyz -cfg config_snapshot.json
```


### Conformational Classification and Ramachandran Analysis

In the present conformational search utilizing the GFN-FF method under vacuum conditions, no distinct local minimum corresponding to the $P_{II}$ (Polyproline II-like) conformation was detected. This computational result is structurally consistent with preceding *ab initio* studies, which indicate that the $P_{II}$ state typically lacks a deep, independent energy basin in the gas phase and is predominantly stabilized by solvation effects (e.g., explicit or implicit water models).

#### Script Specification: `clustering.py`

**Overview**
The `clustering.py` script performs a hybrid structural classification of the conformational minima (nodes) extracted from the `reaction_network.json` file. It maps the geometry of the alanine dipeptide model (Ac-Ala-NHMe) onto a Ramachandran plot and categorizes each structural state into theoretical macrostates using a combination of literature-based heuristic boundaries and a Gaussian Mixture Model (GMM).

**Dependencies**
* Python 3.x
* `numpy`
* `matplotlib`
* `scikit-learn`

**Input Data Requirements**
* `reaction_network.json`: Must contain a `nodes` array with `node_id` and `xyz_file` path parameters.
* Corresponding `.xyz` files containing the optimized atomic coordinates. The script assumes the following atomic index sequence for the backbone definition: C_acetyl [1], N_amide [6], C_alpha [2], C_carbonyl [4], N_amide [7] (0-indexed standard).

**Classification Methodology**
1. **Dihedral Angle Calculation**: Calculates the backbone dihedral angles ($\phi, \psi$) for each optimized geometry. Angles are algorithmically wrapped to the periodic boundary of $[-180^\circ, 180^\circ]$.
2. **Heuristic Region Assignment**: Nodes are deterministically assigned to canonical macrostates ($C_5$, $P_{II}$, $C_7^{eq}$, $\alpha_R$, $\alpha_L$, $C_7^{ax}$) based on bounded rectangular regions derived from structural biology literature. To resolve spatial overlaps, the evaluation priority is strictly ordered (e.g., the $P_{II}$ region is evaluated prior to $C_7^{eq}$).
3. **GMM for Outliers**: Nodes falling outside the predefined heuristic boundaries are temporarily labeled as "Other". A Gaussian Mixture Model with diagonal covariance is applied exclusively to this subset to objectively identify local distributions. The optimal number of mixture components is determined automatically via the Bayesian Information Criterion (BIC).

**Usage**
Execute the script from the directory containing `reaction_network.json` and the corresponding `.xyz` files:
```bash
python clustering.py
```
**Output**
`sramachandran_classified.png`: A 2D scatter plot visualizing the $(\phi, \psi)$ distribution, heuristic boundaries (represented as colored rectangles), and node ID annotations.

Standard Output: A console summary detailing the optimal GMM component count, the corresponding BIC value, and the quantitative distribution of nodes assigned to each macroscopic state.

<img width="3000" height="3000" alt="ramachandran_classified" src="https://github.com/user-attachments/assets/04b469ce-8d87-4ba0-b44d-03bcc1f0aba2" />

`activation_energy_analysis.py`
```
NOTE: The reaction graph contains 2 connected components.
  Component 0: nodes [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]  (main)
  Component 1: nodes [6, 7]  (isolated — excluded)

NOTE: Node(s) [6] assigned to 'C7eq' are isolated and will be excluded from barrier calculations.

Origin       | Destination  |   Barrier (kcal/mol) | Best path (node IDs)
-------------------------------------------------------------------------
C5           | C7eq         |                 1.16 | 0 -> 13 -> 14
C7eq         | C5           |                 1.89 | 14 -> 13 -> 0
C5           | C7ax         |                 2.26 | 0 -> 1 -> 3 -> 4 -> 11
alpha_L      | alpha_R      |                 2.34 | 9 -> 8 -> 20
C7ax         | C5           |                 2.62 | 11 -> 4 -> 3 -> 1 -> 0
C7ax         | C7eq         |                 2.62 | 11 -> 4 -> 3 -> 1 -> 0 -> 13 -> 14
C7eq         | C7ax         |                 2.99 | 14 -> 13 -> 0 -> 1 -> 3 -> 4 -> 11
alpha_R      | alpha_L      |                 3.35 | 20 -> 8 -> 9
alpha_L      | C5           |                 3.78 | 9 -> 8 -> 20 -> 21
alpha_L      | C7eq         |                 4.20 | 9 -> 15 -> 23 -> 27 -> 0 -> 13 -> 14
alpha_L      | C7ax         |                 4.20 | 9 -> 15 -> 23 -> 27 -> 0 -> 1 -> 3 -> 4 -> 11
alpha_R      | C5           |                 4.79 | 20 -> 21
alpha_R      | C7eq         |                 5.21 | 20 -> 8 -> 9 -> 15 -> 23 -> 27 -> 0 -> 13 -> 14
alpha_R      | C7ax         |                 5.21 | 20 -> 8 -> 9 -> 15 -> 23 -> 27 -> 0 -> 1 -> 3 -> 4 -> 11
C5           | alpha_R      |                 7.62 | 21 -> 20
C5           | alpha_L      |                 7.62 | 21 -> 20 -> 8 -> 9
C7ax         | alpha_R      |                 8.41 | 11 -> 4 -> 3 -> 1 -> 0 -> 27 -> 23 -> 15 -> 9 -> 8 -> 20
C7ax         | alpha_L      |                 8.41 | 11 -> 4 -> 3 -> 1 -> 0 -> 27 -> 23 -> 15 -> 9
C7eq         | alpha_R      |                 8.78 | 14 -> 13 -> 0 -> 27 -> 23 -> 15 -> 9 -> 8 -> 20
C7eq         | alpha_L      |                 8.78 | 14 -> 13 -> 0 -> 27 -> 23 -> 15 -> 9

=== Lowest escape barrier per conformational basin ===
  C5          : 1.16 kcal/mol  ->  C7eq
  C7eq        : 1.89 kcal/mol  ->  C5
  alpha_R     : 3.35 kcal/mol  ->  alpha_L
  alpha_L     : 2.34 kcal/mol  ->  alpha_R
  C7ax        : 2.62 kcal/mol  ->  C5

NOTE: All barriers are GFN-FF potential energy barriers (kcal/mol).
      They do not include ZPE, entropy, or thermal corrections and
      must not be interpreted as free-energy barriers.
```

### Kinetic Network Analysis and Basin-to-Basin Activation Barriers

An effective activation barrier analysis was performed to evaluate the kinetic connectivity between the identified conformational basins. The basin-to-basin transition barriers were derived by finding the minimax path (the pathway minimizing the highest transition state energy) between sets of nodes assigned to each conformation basin on the Minimum Spanning Tree (MST) of the reaction network.

#### Key Observations

- **Fast Equilibrium in Extended States:**
  The lowest activation barriers were observed between the $C_5$ and $C_{7\mathrm{eq}}$ states
  (1.16 kcal/mol for $C_5 \to C_{7\mathrm{eq}}$, and 1.89 kcal/mol for $C_{7\mathrm{eq}} \to C_5$),
  suggesting that these extended structures interconvert rapidly and together constitute the dominant
  potential energy minimum of the system under vacuum conditions.


- **Rate-Limiting Steps:**
  Transitions from the extended states ($C_5$, $C_{7\mathrm{eq}}$) to the helical states ($\alpha_R$, $\alpha_L$)
  require substantially higher activation energies (e.g., 7.62 kcal/mol for $C_5 \to \alpha_R$),
  representing the kinetic bottleneck for large-scale conformational reorganization.

- **Potential Energy Preference:**
  The barrier for escaping $\alpha_R$ toward $C_5$ (4.79 kcal/mol) is lower than the reverse (7.62 kcal/mol),
  satisfying detailed balance. This is consistent with the well-established result that extended conformations
  are energetically favored over helical ones in the gas phase owing to intramolecular hydrogen bonding.

#### Computational Limitations

The quantitative barrier values are derived from the GFN-FF force field and represent **potential energy
barriers only**; zero-point energy, entropic contributions, and thermal corrections are not included,
and the values should not be interpreted as free-energy barriers. The accuracy of GFN-FF transition-state
energies has not been systematically validated against high-level *ab initio* references for this system,
and quantitative comparison with experiment should be treated with caution.

A near-zero reverse barrier on one edge (edge 95, $\Delta E_\mathrm{rev} = 0.004$ kcal/mol) was detected;
this is consistent with a nearly barrierless region of the potential energy surface rather than an
optimization artifact, as the stored forward barrier (0.307 kcal/mol) is physically reasonable.
Additionally, isolated nodes (e.g., node 6) lacking TS connections to the main graph component were
automatically excluded from the minimax path calculations.

