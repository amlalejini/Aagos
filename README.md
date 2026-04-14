# Aagos

**A**uto-**A**daptive **G**enetic **O**rganization **S**ystem

[![OSF](https://img.shields.io/badge/data%20%40%20OSF-10.17605%2FOSF.IO%2FWM659-blue)](https://osf.io/wm659/)
[![web-model](https://img.shields.io/badge/go_to-interactive_web--based_model-purple)](https://lalejini.com/Aagos/web/Aagos.html)

## Overview

### Abstract

> Overlapping genes are featured in genetic architectures across the tree of life, yet the prevalence of gene overlap can vary dramatically across species.
  Indeed, there is an evolutionary trade-off in how genetic information is organized.
  Gene overlap allows for more compact information storage and enables tight physical coupling between genes, whereas more modular, non-overlapping arrangements allow constituent genes to be modified independently.
  Such modularity has been shown to facilitate adaptive evolution (i.e., increase evolvability), but how does it come to exist?
  The selective advantage of this evolvability manifests only in the long term, and taken alone provides a weak driver for gene segregation.
  Rapid environmental change has been proposed as a mechanism to heighten the selective importance of evolvability.
  Furthermore, specific types of changing environments have been shown to promote the evolution of functional modularity.
  Here, we expand this theoretical framework to predict that in any rapidly changing environment genetic architectures will be selected that allow independent evolution of organismal traits that respond to independently varying environmental features.
  However, gene segregation produces longer genomes that are larger mutational targets. As such, high mutation rates can produce a counterbalancing selection pressure for more compact genetic architectures.
  We use computational models to verify predictions for how environmental change and mutation rate shape the evolution of gene segregation.
  Specifically, we demonstrate that changing environments promote gene segregation, and we confirm that segregated genes are better able to adapt to novel environments.
  In contrast, we find that high mutation rate is sufficient to countervail the benefits of gene segregation and drive the evolution of gene overlap in order to reduce the mutational load on coding regions.

## Repository guide

- `docs/` - Contains supplemental documentation. E.g., getting started guide for compiling and running the model locally.
- `experiments/` - Contains experiment configuration files, HPC job submission scripts, data analyses, and generated plots.
- `genome-architectures/` - Contains enumerated (meaningfully different) genetic architectures for different gene count, gene size, and genome length parameterizations.
- `hpc-env/` - Contains bash scripts for configuring HPC software environments for running experiments.
- `scripts/` - Contains Python scripts and utilities used for data analyses / managing experiments.
- `source/` - Contains source code for computational model / experiment software.
- `third-party/` - Contains Empirical library dependency as a git submodule.
- `web/` - Contains web build of computational model. Access [here](https://lalejini.com/Aagos/web/Aagos.html).