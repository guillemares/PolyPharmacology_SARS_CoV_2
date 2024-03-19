# PolyPharmacology_SARS-CoV-2

This repository offers a series of Python scripts designed to analyze the three-dimensional structure of multiple RNA molecules.

--- 3D_RNAlysis is currently under development ---

## Features
- **Base Pair Detection:** The included scripts allow for the identification of single chain RNAs intra-bp, returning the type of interaction (Watson-Crick, Hogsteen, Sugar) using Biotite packages.
- **Hydrogen Bond Identification:** In addition to base pairs, the user is allowed to compute all H-bonds in a RNA structure for a more exhaustive analysis.
- **RNA Superimposition:** Given a multiple set of RNAs, 3D-RNAlysis can superpose the structures and compute RMSD and TM-score using US-align packages using two methods.
  - 1 vs 1: Every structure is superimposed with the rest of files in order to generate a matrix of RMSDs and TM-scores
  - All vs All: Given a multiple set of structures, a single superposition is given for all RNAs with a unique average RMSD.
- **Analysis and Clustering:**

## Issues and Suggestions
If you encounter any problems or have any suggestions, please open an issue in this repository. We'll be happy to assist you.


## About the authors
Hi there! 👋 My name is [Guillem Arasa](https://www.linkedin.com/in/guillem-arasa-estivill-1551ab1b3/) and I'm currently pursuing a Master's degree in [Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry](http://www.ub.edu/computational_modelling/) at Universitat de Barcelona (UB) and Universitat Politècnica de Catalunya (UPC). This code is being developed by myself, [Dr. Isaac Filella](https://github.com/IFilella) and [PhD Júlia Vilalta](https://github.com/juliavilmor) at [Barcelona Supercomputing Center](https://www.bsc.es/ca).
