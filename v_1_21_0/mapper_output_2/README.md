This technical note focuses on the operational validation and computational feasibility of the reaction pathway exploration tool using GFN2-xTB. It is not intended to provide a highly accurate quantitative assessment of the chemical reactions.

The results presented here were obtained through exploration with GFN2-xTB.
clustered_network_spatial.png: Cluster 1 mainly consists of the products of the Claisen rearrangement, while Clusters 2 and 3 consist of the reactants.

To rapidly evaluate the operation of the tool, GFN2-xTB was employed. Although semi-empirical methods possess inherent limitations in accuracy for detailed pathway explorations, this approach was chosen specifically to ensure that the test computations could be completed under the hardware constraints.

2026-02-27 21:04:46 [INFO    ] multioptpy.Wrapper.mapper: All candidate (EQ, pair) combinations appear exhausted after 10000 resampling attempts. Stopping.

GFN2-xTB: _J. Chem. Theory Comput._ **2019**, 15, 3, 1652–1671


### Execution Procedure
The reaction pathway exploration was executed using the run_mapper.py script included in the MultiOptPy package. The configuration file (config_snapshot.json) and the input coordinate file (mapper_test.xyz) were placed in the same working directory prior to execution.


#### System Requirements and Prerequisites

The computations were performed under the following software environment:

Operating System: Ubuntu 24.04.3 LTS (Noble Numbat)

Programming Language: Python 3.12

Software Package: MultiOptPy v1.21.0


#### Command:
```
python run_mapper.py mapper_test.xyz -cfg config_snapshot.json
```

