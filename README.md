This repository contains data and analysis files that support "Social brain hypothesis for cooperation in primate social networks," by N. G. MacLaren, L. Meng, M. Collier, and N. Masuda. If you have any difficulty reproducing our results or have questions or input, please don't hesitate to contact the authors or open a new issue in this repository.

## Evaluating the Cooperation Threshold

Allen et al. (2017) provide [MATLAB code](https://zenodo.org/record/276933#.Y9CG6MbMJhE) to calculate the cooperation threshold. We (L. Meng) have implemented their procedure in Python. The file `evaluate-threshold.py` calculates the cooperation threshold for each network in a given directory (in this case, downloaded from the [Animal Social Network Repository](https://github.com/bansallab/asnr)). In the `demo` directory, we provide the same code written as a Python function (`eval_threshold()`) and demonstration of that function's use (`evaluate-threshold-demo.py`) on an example network (`weighted_primate_matrix_5.graphml`). Note that the function can fail on large and/or disconnected networks.

## Dependencies

The file `analysis.R` depends on the [MuMIn](https://CRAN.R-project.org/package=MuMIn) package, available on CRAN. The file `evaluate-threshold.py` depends on [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), and [NetworkX](https://networkx.org/). 
