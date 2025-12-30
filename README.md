# OmniDesign Optimizer (MATLAB)
**Version:** 1.0  
**Status:** Public Release  
**License:** MIT / Academic Free License  
**Principal Investigator:** Prof. Antonio Franchi (University of Twente & Sapienza University of Rome)

---

## Overview
The **OmniDesign Optimizer** is a numerical framework for designing the geometry of fully-actuated and omnidirectional multirotor aerial vehicles. It solves for the optimal orientation of rotors (thrust vectors) to maximize **Kinematic Isotropy**.

The core algorithm uses a **Global Manifold Exhaustion** strategy on the projective space $(\mathbb{RP}^2)^N$. The suite includes advanced topological analysis tools to verify the **Tangent Torus Reduction** and the **N-5 Scaling Law** proposed in the accompanying manuscript.

---

## Video Demonstrations

### Dynamic Manifold Visualization
These animations validate the **N-5 Scaling Law** by simulating the vehicle traversing its solution manifold. They visualize the "Design Nullspace"—a continuous path where the vehicle can reconfigure itself without losing optimality.

**How to read the dashboard:**
* **Top Panel (Physical):** Displays 3 distinct optimal branches ("isomers") morphing in sync alongside random control groups. Notice that while the rotors rotate continuously, their relative phase offsets remain rigorously locked (**Affine Phase Locking**).
* **Bottom Panel (Metrics):** Shows real-time Singular Values ($\sigma_1 \dots \sigma_6$) and Force/Torque Manipulability Ellipsoids.

**Key Observation:**
* **Optimal Branches:** The ellipsoids remain perfectly spherical and static, and singular values stay constant. This proves **Isotropic Invariance**: the vehicle can morph continuously while maintaining maximum control authority.
* **Random/Broken Cases:** The ellipsoids deform ("breathe") and shear, indicating a time-varying loss of control authority.

#### 1. The Cube (N=8)
[PASTE YOUR CUBE VIDEO LINK HERE]

#### 2. The Octagon (N=8)
[PASTE YOUR OCTAGON VIDEO LINK HERE]

---

## Repository Contents

The codebase is organized into three primary modules:

### 1. Core (`src/core/`)
* **`omnidesign_optimizer.m`** (Design Tool)
    * **Usage:** The primary solver. Generates optimal geometries given a chassis shape (Polygon, Platonic Solid, or Custom).
    * **Output:** Saves a `.mat` file containing thousands of converged solutions.
* **`omnidesign_benchmarker.m`** (Validation Tool)
    * **Purpose:** Runs the ablation study comparing Log-Volume ($J_{vol}$) vs. Condition Number ($\kappa$).
    * **Output:** Generates data for scalability tables and plots for the solution landscape (Solution Landscape & Local Minima).

### 2. Analysis (`src/analysis/`)
* **`torus_params_analyzer.m`** (Dimensionality Reduction)
    * **Purpose:** Takes the raw point cloud from the Optimizer and fits a semi-ellipsoid manifold to each rotor's solution set.
    * **Input:** The `.mat` file from the Optimizer.
    * **Output:** Generates `TorusMapping_Data.mat` (Processed angular coordinates).
* **`n5_scaling_law_analyzer.m`** (Topology Verification)
    * **Purpose:** Analyzes the coordination between the intrinsic angles $\theta$ to verify the $K=N-5$ Scaling Law.
    * **Process:** Uses high-dimensional local PCA and clustering to identify discrete isomers (Topological Branches).
    * **Output:** Visualizes disconnected loops and computes the "Cluster Spread" confidence metric.

### 3. Visualization (`src/viz/`)
* **`cube_singular_values_animation.m`**
    * **Demo:** An interactive 3D visualizer for the Cube (N=8) chassis.
    * **Features:** Compares optimal isomers against random configurations and "broken" topologies, displaying real-time manipulability ellipsoids and singular value spectra.
* **`octagon_singular_values_animation.m`**
    * **Demo:** An interactive visualizer for the Planar Octagon (N=8) chassis.
    * **Features:** Demonstrates specific star-polygon winding numbers ($\{8,3\}, \{8,4\}, \{8,5\}$) and their impact on kinematic isotropy.

---

## Workflow for Reproducibility

1.  **Generate Design:** Run `src/core/omnidesign_optimizer.m` for a target chassis (e.g., `config.shape = 'dodecahedron'`).

2.  **Verify Performance:** Run `src/core/omnidesign_benchmarker.m` to reproduce the speed comparison and solution quality validation tables.

3.  **Parametrize:** Run `src/analysis/torus_params_analyzer.m` on the generated optimizer data to map the solution manifold.

4.  **Analyze Topology:** Run `src/analysis/n5_scaling_law_analyzer.m` to confirm the branch count ($K=N-5$) and visualize the solution isomers.

5.  **Visualize:** Run the specific chassis animations in `src/viz/` to inspect the physical meaning of the singular values.

---

## Authors & Affiliations
**  Antonio Franchi**
* [Robotics and Mechatronics (RaM)](https://www.ram.eemcs.utwente.nl/), University of Twente, Enschede, The Netherlands
* [Department of Computer, Control and Management Engineering (DIAG)](https://diag.uniroma1.it/), Sapienza University of Rome, Rome, Italy

---

## Citation
If you use this software in your research, please cite the accompanying publication:

> **Antonio Franchi**, "The N-5 Scaling Law: Topological Dimensionality Reduction in the Optimal Design of Fully-actuated Multirotors," *arXiv preprint arXiv:2512.23619*, 2025.  
> DOI: [10.48550/arXiv.2512.23619](https://doi.org/10.48550/arXiv.2512.23619)

**BibTeX:**
```bibtex
@misc{franchi2025n5,
  title        = {The {N-5} Scaling Law: Topological Dimensionality Reduction in the Optimal Design of Fully-actuated Multirotors},
  author       = {Antonio Franchi},
  year         = {2025},
  eprint       = {2512.23619},
  archivePrefix= {arXiv},
  primaryClass = {cs.RO},
  doi          = {10.48550/arXiv.2512.23619},
  url          = {[https://arxiv.org/abs/2512.23619](https://arxiv.org/abs/2512.23619)}
}
---
