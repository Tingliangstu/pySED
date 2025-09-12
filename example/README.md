# Examples for Using pySED to Calculate the Spectral Energy Density (SED)

This directory contains example cases demonstrating how to use **pySED** for SED calculations on different types of systems.

---

## 📌 Recommended Practice
It is **strongly recommended** to reproduce the results of the provided examples **before** applying pySED to your own system.  
This ensures you are familiar with the workflow and can verify your setup.

---

## 📂 Example Categories

### 1D Systems
- **`CNT`**  
  Example for a one-dimensional system (Carbon Nanotube).  
  Use this to learn how to set up and analyze SED for 1D materials.

### 2D Systems
- **`In_plane_graphene_gpumd`**  
  Example for in-plane graphene using GPUMD.  
  
- **`MoS2_gpumd`**  
  Example for monolayer MoS<sub>2</sub> using GPUMD.  
  
  These are suitable for learning SED analysis of two-dimensional materials.

### 3D Systems
- **`Silicon_primitive_gpumd`**  
  Example for bulk silicon.  
  Recommended for learning SED analysis of three-dimensional crystalline materials.
  
  
### tutorials
- **`tutorials`**  
  Jupyter Tutorial for MoS<sub>2</sub>.  

---

## 🔍 Reference for Method Validation
- **`Ref_Phonon_dispersion_from_phonopy`**  
  Contains phonon dispersion calculated using the **Lattice Dynamics (LD)** method via **phonopy**.  
  This can be used to validate the reliability of the SED method by comparing results.

---

## 💡 Tips
- Start with the **Graphene** or **Silicon** examples to get familiar with the workflow.
- Use the reference phonon dispersion to cross-check your SED results.
- Adjust parameters carefully when moving from example systems to your own material.

---

**Author:** pySED Development Team (Ting Liang and Wenwu Jiang) 

**Version:** v2.2.0 and above