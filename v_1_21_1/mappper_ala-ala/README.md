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

