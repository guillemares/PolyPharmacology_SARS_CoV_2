# 3D-RNAlysis

This repository offers a series of Python scripts designed to analyze the three-dimensional structure of multiple single strain RNA molecules.

------ **3D_RNAlysis is currently under development** ------

## Features
- **Base Pair Detection:** The included scripts allow for the identification of single chain RNAs intra-bp, returning the type of interaction using Biotite packages.
  - *Glycosidic bond:* Cis - Trans.
  - *Edge type:* Watson-Crick - Hoogsteen - Sugar - Wooble GU.
  
- **Hydrogen Bond Identification:** In addition to base pairs, the user is allowed to compute not only base pair H-bonds but all H-bonds in a RNA structure for a more exhaustive analysis.

- **RNA Superimposition:** Given a multiple set of RNAs, 3D-RNAlysis can superpose the structures and compute RMSD and TM-score using US-align packages using two methods.
  - *1 vs 1:* Every structure is superimposed with the rest of files in order to generate a matrix of RMSDs and TM-scores
  - *All vs All:* Given a multiple set of structures, a single superposition is given for all RNAs with a unique average RMSD.

- **Analysis and Clustering:** Using the information calculated in the RNA Superimposition, 3D-RNAlysis provides two possible classifications.
  - *Heatmaps*
  - *Clustering*

## Prerequisites
- Linux with Ubuntu version (20.04.5 LTS or higher).
- Python3 (3.8.17) with sys, numpy, pandas, seaborn, matplotlib, glob, os, scipy packages installed
- Biotite (0.39.0) ---> See [installation guide](https://www.biotite-python.org/install.html).
- USAlign (Version 20231222) ---> See [installation guide](https://zhanggroup.org/US-align/).

## Installation
## Usage

## Issues and Suggestions
If you encounter any problems or have any suggestions, please open an issue in this repository. We'll be happy to assist you.


## About the authors
Hi there! 👋 My name is [Guillem Arasa](https://github.com/guillemares/) and I'm currently pursuing a Master's degree in [Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry](http://www.ub.edu/computational_modelling/) at Universitat de Barcelona (UB) and Universitat Politècnica de Catalunya (UPC). This code is being developed by myself, [Dr. Isaac Filella](https://github.com/IFilella) and [PhD Júlia Vilalta](https://github.com/juliavilmor) at [Barcelona Supercomputing Center](https://www.bsc.es/ca).

[![Linkedin](https://img.shields.io/badge/-LinkedIn-blue?style=flat-square&logo=Linkedin&logoColor=white&link=enlace_al_perfil_de_LinkedIn)](https://www.linkedin.com/in/guillem-arasa-estivill-1551ab1b3/)
[![GitHub](https://img.shields.io/badge/-GitHub-181717?style=flat-square&logo=github&logoColor=white&link=enlace_al_perfil_de_GitHub)]([enlace_al_perfil_de_GitHub](https://github.com/guillemares/)https://github.com/guillemares/)
