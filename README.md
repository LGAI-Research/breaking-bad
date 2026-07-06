# Breaking Bad: Component-Wise Parent Deletion for Score-Based Causal Discovery

This repository provides the official implementation for our paper:

> *Breaking Bad: Component-Wise Parent Deletion for Score-Based Causal Discovery*

We propose **parent deletion**, a novel perturbation operator for score-based causal discovery that delete all incoming edges from the parents of a target node (or an entire undirected component) simultaneously. 
This operator is theoretically sound, computationally efficient, and seamlessly integrates with existing score-based methods via an Iterated Local Search (ILS) framework. Extensive experiments on synthetic and realistic datasets demonstrate consistent improvements across a wide range of settings.
[[paper]](https://openreview.net/forum?id=oxdOxCxgCu)

---

## Public Release Notice
This public release omits several files that directly depend on external libraries due to licensing and intellectual property considerations.

Consequently, experiments cannot be executed directly using the repository as-is. To reproduce the full experimental pipeline, users need to obtain the corresponding open-source implementations separately and place them in the appropriate locations.

Detailed references to all related open-source projects are provided in the **Reference** section and main paper, including the repositories from which the required functionality can be obtained.

We gratefully acknowledge the contributions of the original authors.

---
# Repository Overview

This repo provides:

- C++ implementation of score-based search
- Support for various baselines (GES, OPS, XGES, LGES, BOSS, PC)
- Parent-deletion perturbation operators
- Experimental pipeline (ER / Scale-Free graphs and realistic experiments)
- Automatic evaluation and summary generation
- Multi-seed aggregation (default: 30 runs)

---

## Table of Contents

- [Setup](#setup)
  - [Python Environment](#python-environment)
  - [Building the C++ Search Engine](#building-the-c-search-engine)
- [Data Generation](#data-generation)
  - [Erdős–Rényi (ER) Graphs](#erdősrényi-er-graphs)
  - [Scale-Free (SF) Graphs](#scale-free-sf-graphs)
  - [Biologically motivated benchmark dataset](#biologically-motivated-benchmark-dataset)
- [Running Experiments](#running-experiments)
  - [Full Experiment Pipeline](#full-experiment-pipeline)
  - [Single Experiment Example](#single-experiment-example)
  - [Configuration](#configuration)
- [Method Variants](#method-variants)
- [Experiment Pipeline Overview](#experiment-pipeline-overview)
- [Reference](#reference)

---

# Using the Python Package

Create a conda environment:

```bash
conda create -n <env_name> python=3.10
conda activate <env_name>
pip install -r requirements.txt
```

Note: If you encounter issues with rpy2, Installing rpy2 may require R to be installed and accessible in your system PATH.

---

# Building the C++ Code

Install dependencies:

```bash
conda install -c conda-forge boost openblas r-bnlearn
```

Clone required libraries:

```bash
cd breaking-bad
mkdir -p lib
cd lib
git clone https://github.com/gabime/spdlog.git
git clone https://github.com/rogersce/cnpy.git
```

## BIC scorer

Graphs are scored with a Gaussian BIC scorer whose sources (`BICScorer.h` / `BICScorer.cpp`) are obtained separately from open-source projects rather than shipped here. To build:

1. Place those two files into `breaking-bad/src/`.
2. Complete `breaking-bad/src/BICScoreloader.cpp` following the short notes at the top of that file.
3. The scorer needs a matrix library (such as Eigen). Install it yourself and make it available to the build: add its `find_package()`/link in `breaking-bad/src/CMakeLists.txt` where the comment indicates, or pass its include path at configure time (e.g. `cmake -DCMAKE_CXX_FLAGS="-I<path-to-matrix-library-include>" ..`).

`main.cpp` reaches the scorer only through a small interface, so no other file needs changes.

Build:

```bash
cd breaking-bad
mkdir build
cd build
cmake ..
make
```

Note: If configuration fails in `lib/cnpy` due to “CMake < 3.5”, rerun CMake with `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`.

---

# Data Generation

We consider both synthetic and realistic benchmark datasets.

For synthetic experiments, we generate datasets under two graph families:

- **ER (Erdős–Rényi)**
- **SF (Scale-Free)**

For realistic evaluation, we use a biologically motivated benchmark network.

---

## ER Graphs

Script:
```
paper_experiments/simple_exp/simple_simulation.py
```

- DAGs are generated using `dag_avg_deg(d, ρ, ...)` from the `sempler` library.
- Data is sampled from linear-Gaussian SEMs: `v = w_V^T pa_V + ε_V`.
- Structure controlled by `graph_seed`; data controlled by a seed derived from `graph_seed + n_samples`.

Output:


| Type | Path |
|------|------|
| Ground truth (CSV) | `data/ground_truth/True_DAG_var<d>_avg_deg<ρ>_graph_num<seed>.csv` |
| Samples (NPY) | `data/samples/Data_var<d>_avg_deg<ρ>_n_samples<n>_graph_num<seed>.npy` |


---

## Scale-Free (SF) Graphs

Script:
```
paper_experiments/simple_exp/create_scale_free.py
```

- Undirected graphs generated via `nx.barabasi_albert_graph(d, m)`.
- A random node ordering is used to orient edges, ensuring a valid DAG.
- Data is sampled from the same linear-Gaussian SEM.


Output:

| Type | Path |
|------|------|
| Ground truth (CSV) | `data/ground_truth/ScaleFree_DAG_var<d>_seed<seed>.csv` |
| Samples (NPY) | `data/samples/ScaleFree_Data_var<d>_n<n>_seed<seed>.npy` |

---


## Biologically motivated benchmark dataset

Script:
```
paper_experiments/simple_exp/create_realworld.py
```

- A biologically motivated benchmark network originally constructed from high-dimensional plant gene expression data.
- Leveraged from the GeneNet package and distributed through the bnlearn repository.
- The true network structure is loaded from a local rds file.
- Data generation is performed using the bnlearn R package (via rpy2) to ensure fair evaluation according to the original network file configurations.



Output:

| Type | Path |
|------|------|
| Ground truth (CSV) | `True_DAG_ARTH150.csv` |
| Samples (NPY) | `data/Realworld/ARTH150/Data_ARTH150_<n>_<seed>.npy` |

---


# Running Large-Scale Experiments

All experiments are controlled by:

```
paper_experiments/run_all_experiments.sh
```

Run:

```bash
cd paper_experiments
bash run_all_experiments.sh
```

---

# Single Experiment Example

```bash
cd paper_experiments
bash exp_run.sh --var 20 --deg 2.0 --samples 1000 --alpha 2.0
```

---

## Experiment Configuration

Edit `run_all_experiments.sh` to configure the experiment grid. The key parameters are:

| Parameter | Description | Example values |
|-----------|-------------|----------------|
| `VAR_LIST` | Number of nodes (*d*) | `(50 100 200)` |
| `DEG_LIST` | Average degree (*ρ*) for ER graphs | `(2.0 3.0 5.0)` |
| `SAMPLES_LIST` | Number of samples (*n*) | `(1000 10000 100000)` |
| `ALPHA_LIST` | BIC penalty hyperparameter (*λ*) | `(1 2 4)` |
| `VARIANTS_LIST` | Method variant ID (see [Method Variants](#method-variants)) | `(0 1 2 3)` |
| `BASELINE_LIST` | Base search algorithm | `("ges" "ops" "xges" "boss" "lges-safe" "lges-cons")` |
| `GRAPH_TYPE` | Graph family | `"ER"` or `"SF"` |
| `RUNS` | Number of random seeds | `30` |
| `M_LIST` | Scale-free density parameter (*m*) | `(2 4)` |

**Results** are saved under:

```
paper_experiments/outputs/<baseline>/<graph_type>/var<d>/deg<ρ>/
```

**Summary statistics** (mean ± std over 30 seeds) are saved under:

```
paper_experiments/outputs/summary/<baseline>/<graph_type>/var<d>/deg<ρ>/
```

---


## Method Variants

The `--variant` flag controls the perturbation strategy applied on top of each baseline:

| Variant ID | Name suffix | Perturbation Strategy |
|:----------:|-------------|----------------------|
| 0 | *(base)* | No ILS perturbation (vanilla baseline) |
| 1 | `-D` | Single-edge deletion perturbation (XGES-style) |
| 2 | `-DP` | **Component-wise parent deletion** (Alg. 2 in the paper) |
| 3 | `-DP+` | **Three-phase scheduled perturbation** (DP → single-node → single-edge) |

For example, `--baseline ges --variant 2` runs **GES-DP**, and `--baseline boss --variant 3` runs **BOSS-DP+**.

---

## Experiment Pipeline

```
run_all_experiments.sh            ← Grid search over all configurations
    ↓
batch_run.py                      ← Multi-seed runner, result aggregation
    ↓
exp_run.sh                        ← Single-seed experiment launcher
    ↓
build/src/Search (C++ binary)     ← Core search algorithm execution
    ↓
simple_exp.py                     ← Evaluation (SHD, F1, ecSHD, etc.)
```

---


# Running Real-World Experiments

Real-world experiments (e.g., ARTH150) use a dedicated launcher.
To run real-world experiments:

```bash
cd paper_experiments
bash exp_run_realworld.sh
```

This script runs the search algorithm on the biologically motivated benchmark dataset and performs evaluation in the same manner as synthetic experiments.

---


## BOSS Baseline Execution

Runs that use BOSS initialization (`--baseline boss`) require a BOSS-produced
graph as the starting structure. This step is **not included** in the public
release: BOSS is provided by the open-source **Tetrad** library, which is not
distributed here.


You must therefore generate the initial graph yourself using an open-source BOSS
implementation **before** launching those experiments


---

# Reference

This project builds upon and reuses components from several open-source repositories.
We thank the original authors for making their implementations publicly available.

## Core Structure Learning Implementations

- **XGES**: https://github.com/ANazaret/XGES  
The implementation of score computation is adapted from:
  - https://github.com/ANazaret/XGES/blob/main/src-cpp/BICScorer.cpp
  - https://github.com/ANazaret/XGES/blob/main/src-cpp/BICScorer.h

- **LGES**: https://github.com/CausalAILab/lges  

- **Tetrad / Py-Tetrad**: https://github.com/cmu-phil/tetrad / https://github.com/cmu-phil/py-tetrad  


## Python Libraries

- **causal-learn**: https://github.com/py-why/causal-learn  

- **sempler**: https://github.com/juangamella/sempler  

- **gadjid**: https://github.com/CausalDisco/gadjid  

## Network Datasets

- **GeneNet (R package)**: https://github.com/cran/GeneNet  
  
- **bnlearn**: https://www.bnlearn.com/bnrepository/  

---

## Citation

If you use this code, please cite the following paper:

```
@inproceedings{park2026breaking,
  title={Breaking Bad: Component-Wise Parent Deletion for Score-Based Causal Discovery},
  author={Park, Min Woo and Yun, Taehui and Jang, YoungIn and Yeom, Yoonseok and Kim, Jonghwan and Kang, Jiyeon and Kim, Songseong and Jung, Hyemin and Lee, Sangmin and Jang, Jongseong and Lee, Sanghack},
  booktitle={Proceedings of the Forty-Second Conference on Uncertainty in Artificial Intelligence},
  year={2026},
}
```