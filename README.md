# Adaptive TT-Rounding using KRP Structure

This repository contains MATLAB code to reproduce the experiments presented in the paper:

**"Adaptive Randomized Tensor Train Rounding using Khatri-Rao Products"** by *Hussam Al Daas, Grey Ballard, Laura Grigori, Mariana Martinez Aguilar, Arvind K. Saibaba  and Bhisham Dev Verma* (Submitted, 2025).

---

## 📁 Repository Structure & Usage

| File/Folder             | Description |
|-------------------------|-------------|
| `NormEstimation.m`      | Reproduces the norm estimation results (Figure 3, Section 3.4). |
| `SyntheticTest.m`       | Runs the synthetic experiments (Section 4.1). Depends on `TTcore/` and `TTrandomized/`. |
| `MaternTest.m`          | Executes the Matérn kernel-based experiments (Section 4.2). Requires `TTcore/`, `TTrandomized/`, and `Dataset/`. |
| `TTGMRESTest.m`         | Runs the GMRES experiments (Section 4.3). Uses `readb.m`, `readcoomat.m`, and `timed_TTGMRES.m`. Depends on `TTcore/`, `TTrandomized/`, and `Dataset/`. |
| `Experimental_results/` | Contains saved figures and data corresponding to experimental results in the paper. |
| `tensor_toolbox/`       | Bundled version of **Tensor Toolbox 3.6** (no external installation required). |

---

