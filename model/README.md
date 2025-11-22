## FEKO models in this directory

This directory contains FEKO models. Files with the `.cfx` extension can be opened and edited directly in FEKO.

- `Horn_Liang_Script.cfx`  
  A horn antenna with **five propagating modes** in the frequency range **3.2 GHz – 3.8 GHz**.

- `dipole_FEM_MoM_hyb.cfx`  
  A dipole antenna **embedded in a cylindrical dielectric cavity**, modelled using a hybrid **FEM–MoM** formulation in FEKO.

> **Note:**  
> These `.cfx` files are **not required** when running the MATLAB scripts.  
> They are mainly used to generate the corresponding `.cfm` and `.pre` files.  
> Additionally, the `.cfx` file determines **which numerical method** (MoM, FEM, MLFMM, hybrid FEM–MoM, etc.) is used in each region of the structure.

---

## Required FEKO files

For each antenna model, FEKO should export:

- a `.cfm` file (mesh/configuration file)
- a `.pre` file (pre-processing / model definition)

These two files **must reside in the same directory** for the MATLAB–FEKO co-simulation scripts to function correctly.

---

## Extending to custom antenna types

To extend this framework to your own antenna geometries, please follow these guidelines when setting up the FEKO model:

1. **Use waveguide ports**  
   - Antenna ports must be defined as **wave ports** in FEKO.  
   - (At the moment, **lumped ports are not supported**.)

2. **Assign a waveguide source to each wave port**  
   - Each wave port must have a **waveguide source** assigned.  
   - For **multimode** structures, it is sufficient to excite **one mode** at each port;  
     the scripts will automatically determine how many propagating modes exist over the frequency band.

3. **Multiple ports are supported**  
   - Antennas with **multiple ports** are supported.  
   - Each port must:
     - be defined as a wave port, and  
     - have a waveguide source assigned.

4. **Avoid unnecessary excitations and monitors**  
   - Do **not** define any additional excitations (plane waves, voltage sources, etc.) or port/field monitors.  
   - Only geometry, materials, wave ports, and their waveguide sources should remain.

---

## Contact

For questions, comments, or bug reports, please contact:

**Chenbo Shi**  
University of Electronic Science and Technology of China (UESTC)  
📧 **chenbo_shi@163.com**