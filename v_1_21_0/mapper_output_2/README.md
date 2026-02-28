This technical note focuses on the operational validation and computational feasibility of the reaction pathway exploration tool using GFN2-xTB. It is not intended to provide a highly accurate quantitative assessment of the chemical reactions.


### Case Study: Claisen Rearrangement
The validation was performed on a Claisen rearrangement system. The resulting network was analyzed using unsupervised clustering based on activation barriers:

Result Summary: The exploration successfully identified the core reaction network.

**clustered_network_spatial.png**: Cluster 1 mainly consists of the products of the Claisen rearrangement, while Clusters 2 and 3 consist of the reactants.


### Execution Procedure
#### 1. System Requirements and Prerequisites
The computations were performed under the following software and hardware environment. It should be noted that the Linux environment was operated under a Microsoft hypervisor (e.g., WSL or Hyper-V).

- Operating System: Ubuntu 24.04.3 LTS (Noble Numbat)
- Hardware Specifications:
  - CPU: AMD Ryzen 7 5800H (8 cores, 16 threads)
  - Memory (RAM): 8 GB (Available: approx. 7.8 GB)
  - GPU: Not utilized (CPU-only execution)
- Programming Language: Python 3.12
- Software Package: MultiOptPy v1.21.0

#### 2. Execution Procedure
To execute the reaction pathway exploration, the configuration file (config_snapshot.json) and the input coordinate file (mapper_test.xyz) were placed within the same working directory. The tool was invoked using the run_mapper.py script provided by the MultiOptPy package.

##### Command:
```
python run_mapper.py mapper_test.xyz -cfg config_snapshot.json
```
#### 3. Computational Performance

To rapidly evaluate the operation of the tool, GFN2-xTB was employed. Although semi-empirical methods possess inherent limitations in accuracy for detailed pathway explorations, this approach was chosen specifically to ensure that the test computations could be completed under the hardware constraints.

Under the specified hardware constraints and utilizing the GFN2-xTB method, the wall-clock time required to complete the pathway exploration test was approximately 72 hours. This metric is provided solely to demonstrate the operational feasibility of the workflow on a standard laptop environment.


Completion Status: The process reached an exhaustive state after 10,000 resampling attempts, successfully mapping the stable nodes and transition states.
```
2026-02-27 21:04:46 [INFO    ] multioptpy.Wrapper.mapper: All candidate (EQ, pair) combinations appear exhausted after 10000 resampling attempts. Stopping.
```
### Identification of Transition State Conformers
- During the exploration of the Claisen rearrangement, the mapping algorithm identified four transition state geometries. These states correspond to chair and boat topologies, including their respective enantiomeric pairs:

- Chair-like Transition States: A chair-like geometry (TS000000.xyz) and its enantiomer (TS000061.xyz) were located. The activation barrier from the reactant to the product for this pathway is $ΔE^{‡}$ = 19.78 kcal/mol at the GFN2-xTB level. This result is close from a numerical perspective to the general model where the chair conformation provides a lower-energy pathway.

- Boat-like Transition States: A boat-like transition state (TS000022.xyz) and its enantiomer (TS000074.xyz) were identified. The corresponding activation barrier for this pathway was calculated as $ΔE^{‡}$ = 27.26 kcal/mol. This relative increase structurally corresponds to the higher-energy nature typically associated with boat conformations.

### Conclusion and Limitations
The reaction pathway exploration using GFN2-xTB reproduced the qualitative trend wherein the chair-like pathway exhibits a lower activation barrier compared to the boat-like pathway. 

However, the interpretation of the numerical results requires strict limitations. The calculated electronic activation barrier ($\Delta E^{\ddagger}$) for the primary chair-like transition state is 19.78 kcal/mol. This value is numerically lower than the experimental macroscopic activation free energy ($\Delta G^{\ddagger}$ ≈ 33.3 kcal/mol at 469.1 K) reported for the Claisen rearrangement of allyl vinyl ether (_J. Am. Chem. Soc._ 1950, 72, 7, 3155–3159). This discrepancy reflects the inherent accuracy limitations of semi-empirical approaches and the absence of thermodynamic corrections (e.g., zero-point energy and entropy contributions) in the present raw $\Delta E^{\ddagger}$ values.


### References
[1] C. Bannwarth, S. Ehlert, and S. Grimme, _J. Chem. Theory Comput._ 2019, 15, 3, 1652–1671.
[2] F. W. Schuler and G. W. Murphy, _J. Am. Chem. Soc._ 1950, 72, 7, 3155–3159.




