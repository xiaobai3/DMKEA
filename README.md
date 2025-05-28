# DMKEA
Dynamic Multi-Knowledge Evolutionary Algorithm for Sparse Large-Scale Multi-Objective Optimization

Author: Lidan Bai (2024)
Based on PlatEMO Framework: https://github.com/BIMK/PlatEMO
- The directory `Problems/` refers to PlatEMO's internal structure and is not included in this repository.
- Benchmark problems such as `Sparse_NN` and `Sparse_SR` used in our experiments are implemented within PlatEMO under `Problems/Multi-objective optimization/`, and you may need to integrate them separately if not already present.

Description
-----------
DMKEA is a sparse multi-objective evolutionary algorithm designed for solving large-scale problems where only a small subset of variables are relevant.

The algorithm introduces:
- Three guidance vectors to estimate variable importance:
  - pv (prior vector) — learned via multi-interval sampling at initialization.
  - fv (filter vector) — refined from top-performing individuals.
  - sv (statistical vector) — continuously updated during evolution via a moving average vote.
- Two-stage operator strategy, both involving binary and real-valued optimization:
  - Stage I: Uses pv and fv to guide mask generation and applies adaptive real-coded variation on active variables.
  - Stage II: Uses sv to guide mask evolution and applies real-coded variation with a fixed mutation rate on active variables.
An  adaptive switching mechanism controls the transition between stages based on search progress, effectively balancing exploration and exploitation.

Structure
---------
- `DMKEA.m` : Entry point for running the DMKEA algorithm.
- `Operator_PvFv.m` : Binary operator guided by pv and fv + real-valued operator with dynamic mutation rate.
- `Operator_Sv.m` : Binary operator guided by sv + real-valued operator with fixed mutation rate.
- `RealVariation_MutateRateByActive.m` : Real-coded variation operator with dynamic mutation rate.
- `RealVariation_MutateRateFixed.m` : Real-coded variation operator with fixed mutation rate.
- `SPEA2_EnvironmentalSelection.m` : Environmental selection used by DMKEA.
- `Problems/` : (Not included here) Benchmark and real-world sparse problems provided by the PlatEMO framework.

Requirements
------------
- MATLAB R2020b or later
- PlatEMO framework (already assumed integrated)

Usage
-----
1. Clone or download this project into your `PlatEMO/Algorithms/Multi-objective optimization/` directory.
2. From MATLAB, run:

   ```matlab
   main_DMKEA
   
Notes
----
1. Problem settings, sparsity levels, and dimensionality can be configured directly in the main_DMKEA.m script.
2. Real-world problem instances may require updating the pv vector during evolution, unlike benchmark problems which use a fixed pv computed at initialization.

Citation
--------
If you use DMKEA in your research, please cite:

L. Bai et al., "Dynamic Knowledge-Guided Evolutionary Algorithm for Sparse Large-Scale Multi-Objective Optimization," 2025. (Under review)

Acknowledgement
----------------
This algorithm was implemented and tested using the PlatEMO platform:

Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, "PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective Optimization [Educational Forum]," IEEE Computational Intelligence Magazine, vol. 12, no. 4, pp. 73–87, 2017.

License
-------
This code is provided for academic and research use only. Please retain this notice in all copies or redistributions of the code.


