# Generalized Scattering Matrix (GSM) Framework  
## Modeling Implantable Antennas in Multi-layered Spherical Media

---

## Table of Contents
1. [Overview](#overview)
2. [Software Requirements](#software-requirements)
3. [Precomputation Scripts](#precomputation-scripts)
4. [Example Scripts](#example-scripts)
5. [FEKO Models](#feko-models)
6. [Requirements for Extending to Custom Antennas](#requirements-for-extending-to-custom-antennas)
7. [Utilities](#utilities)
8. [Citation](#citation)
9. [Contact](#contact)

---

## Overview

This repository provides MATLAB + FEKO implementations for evaluating the electromagnetic performance of antennas embedded in **multi-layered spherical media**.  
All examples are based on the **GSM–SSO (Generalized Scattering Matrix + Spherical Scattering Operator)** framework proposed in [1].

The repository includes complete SSO computations and GSM-coupling procedures used in the main examples of [1].

> ⚠️ **Important Notice**  
> These codes are intended to demonstrate and reproduce the procedures described in the paper.  
> They are **not** general-purpose EM tools.  
> Please use, modify, and extend them cautiously and at your own risk.

This repository is released as supplementary material for [1].  
You may cite [1] when referencing this repository.

---

## Software Requirements

### ✔ FEKO
Tested version:
- **Altair FEKO 2021.1-5697 (x64)**

The scripts will likely function in later FEKO versions, but these have not been fully validated.

> **Environment Variable Requirement**  
> Add FEKO’s `bin` folder to the system PATH:
> ```
> C:\Program Files\Altair\2021.1\feko\bin
> ```

### ✔ MATLAB
- **MATLAB R2024a** (verified)

---

## Precomputation Scripts  
### (These may need to be executed before running the main examples)

These scripts generate the **free-space generalized scattering matrix (GSM)** of the antenna elements.  
Precomputed results are already stored in the `preserved_data/` directory.

---

### 1) Five-mode horn antenna (MoM)
`Co_FEKO_MATLAB_Ex0_Horn_antenna.m`

- Uses FEKO's **MoM solver**  
- Computes the GSM of a horn antenna with **five propagating modes** in **3.2–3.8 GHz**

---

### 2) Dipole encapsulated in a dielectric cylindrical cavity (FEM–MoM)
`Co_FEKO_MATLAB_Ex1_dipole_FEM_MoM_hyb.m`

- Uses FEKO's **FEM–MoM hybrid solver**  
- Computes the GSM of a dipole inside a dielectric cylindrical cavity  

All precomputed data are stored under:

preserved_data/
 ├── CoSim_Ex0_Horn_Liang_Script/
 └── CoSim_Ex1_dipole_FEM_MoM_hyb/

---

## Example Scripts (Directly Runnable)

These scripts correspond to the numerical examples in [1] and can be executed immediately:

| Script                                            | Description                                                  | Paper Section |
| ------------------------------------------------- | ------------------------------------------------------------ | ------------- |
| **Ex1_Homogeneous_Isotropic_Medium.m**            | Horn antenna embedded in a homogeneous isotropic sphere      | Sec. IV-A     |
| **Ex2_Homogeneous_Anisotropic_Medium.m**          | Horn antenna embedded in a homogeneous anisotropic sphere    | Sec. IV-B     |
| **Ex3_Radially_Piecewise_Homogeneous_Medium.m**   | Horn antenna embedded in a radially piecewise homogeneous sphere | Sec. IV-C     |
| **Ex4_Radially_Piecewise_Inhomogeneous_Medium.m** | Horn antenna embedded in a radially piecewise inhomogeneous sphere (ODE-based SSO) | Sec. IV-D     |
| **Ex5_Flexibility_of_the_GSM_SSO_Framework.m**    | Dipole in a cylindrical dielectric cavity embedded in spherical media | Sec. V-C      |

All these examples use GSM data from `preserved_data/`.

---

## FEKO Models

This repository includes FEKO `.cfx` model files located in:

model/

The `.cfx` files:

- define the **geometry**,  
- specify **which numerical method** (MoM, FEM, MLFMM, hybrid, etc.) is used in each region,  
- can be opened and edited directly in FEKO.

### Included models

- **Horn_Liang_Script.cfx**  
  Five-mode horn antenna operating from 3.2 to 3.8 GHz.

- **dipole_FEM_MoM_hyb.cfx**  
  Dipole antenna encapsulated in a cylindrical dielectric cavity (FEM–MoM hybrid).

> These `.cfx` models are **not required when running the MATLAB scripts**.  
> They are only used to generate the `.cfm` and `.pre` files, which **must be kept in the same directory** when running simulations.

---

## Requirements for Extending to Custom Antennas

When using your own FEKO antenna model, ensure:

1. **Each antenna port must be a waveguide port**  
   (lumped ports are not supported at this stage).

2. **Each waveguide port must have a waveguide source**  
   - For multi-mode ports, **only one mode needs to be excited**.  
   - The scripts automatically determine the number of propagating modes across frequency.

3. **Multi-port antennas are supported**  
   Each port must have:
   - a waveguide port, and  
   - an associated waveguide source.

4. **Avoid additional excitations or monitors**  
   - Do **not** include plane waves, lumped sources, or field monitors.  
   - Only geometry, materials, ports, and waveguide sources should remain.

---

## Utilities

The folder:

utilities/

contains functions required for SSO computation, including:

- Analytical (closed-form) SSO evaluation  
- Numerical ODE-based SSO evaluation  
- Spherical wave functions and mode indexing tools  
- Helper utilities used throughout the examples  

---

## Citation

If you use this repository for research, please cite:

**[1]** Shi, Chenbo, et al. "Generalized Scattering Matrix Framework for Modeling Implantable Antennas in Multilayered Spherical Media." *arXiv preprint arXiv:2507.13119* (2025).

---

## Contact

For questions, comments, or collaboration requests:

**Chenbo Shi**  
University of Electronic Science and Technology of China (UESTC)  
📧 **chenbo_shi@163.com**

---