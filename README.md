# OmniDesign Optimizer (MATLAB)
**Version:** 4.1  
**Status:** Public Release  
**License:** MIT / Academic Free License  

---

## Overview
The **OmniDesign Optimizer** is a numerical framework for designing the geometry of fully-actuated and omnidirectional multirotor aerial vehicles. It solves for the optimal orientation of rotors (thrust vectors) to maximize **Kinematic Isotropy**.

The core algorithm uses a **Global Manifold Exhaustion** strategy on the projective space $(\mathbb{RP}^2)^N$. The suite includes advanced topological analysis tools to verify the **Tangent Torus Reduction** and the **N-5 Scaling Law** proposed in the accompanying manuscript.

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

## Citation
If you use this software, please cite the accompanying publication:
> [Authors], "Global Manifold Exhaustion for Omnidirectional Multirotor Design," [Journal/Conference], 2025.
