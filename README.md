[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)
![pySED Logo](https://github.com/Tingliangstu/pySED/blob/main/docs/source/_static/logo.png)

# pySED

**To implement the [phonon SED (spectral energy desity) method in 2010](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), phonon lifetime can be easily calculated.** 

## Main features

- Output a very nice SED image 
- Using Lorentz to fitting the SED peak 
- Fit all peaks at once and output phonon lifetimes, or can fit them individually
- Phonon lifetime can only be obtained qualitatively
- Can be used in parallel and get SED results quickly
- Can interface with [GPUMD](https://github.com/brucefan1983/GPUMD?tab=readme-ov-file) and [LAMMPS](https://www.lammps.org/#gsc.tab=0), capable of capturing quantum dynamics.

## Installation instructions

1)  Download by git or from [source code](https://github.com/Tingliangstu/pySED). **Installing from source is highly recommended, as it is frequently maintained**.
```bash
git clone https://github.com/Tingliangstu/pySED.git
```
2) Install
```bash
cd pySED
pip install .
```
  or one can use:
```bash
python setup.py install --user --prefix=
```

3) If one prefers a one-line installation without cloning manually:
```bash
pip install git+https://github.com/Tingliangstu/pySED.git
```

For convenience, you may want to copy (or link) the files inside scripts folder to a location included in $PATH environment variable.

4) Verify installation
```bash
pySED -h
```
  or
```bash
pysed -h
```
This means that `pySED` or `pysed` can be run directly from the command line.
You can always use the `-h` flag to explore available options and understand how to prepare input files for `pySED`.

## Citations

| Reference             | Cite for what?                    |
| --------------------- | --------------------------------- |
| [1]                   | for any work that used `pySED`    |
| [2]                   | fundamental theory on phonon SED |

## References

[1] Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, Jianbin Xu, [PYSED: A tool for extracting kinetic-energy-weighted phonon dispersion and lifetime from molecular dynamics simulations](https://doi.org/10.1063/5.0278798). J. Appl. Phys. **138**, 075101 (2025).

[2] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, A. J. H. McGaughey, [Predicting phonon dispersion relations and lifetimes from the spectral energy density](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), Phys. Rev. B **81**, 081411 (2010).
	
## Publications using pySED

Please visit the [publications page](https://github.com/Tingliangstu/pySED/tree/main/publications).


## Usage

The [example files](https://github.com/Tingliangstu/pySED/tree/main/example) have so many case provided to show how to use **pySED**. 
If you can successfully reproduce these cases, you are ready to apply `pySED` to your own research systems.  

Users can always use the `pysed -h` flag to explore available options and understand how to prepare input files for pySED.

We also strongly encourage readers to browse through the [pySED paper](https://doi.org/10.1063/5.0278798), 
as it contains a wealth of details that can greatly aid in understanding both the theoretical foundations 
and the practical aspects of the SED methods it presents.


**Online manual** [https://pysed.readthedocs.io](https://pysed.readthedocs.io/en/latest/)


### 1D Systems
- **Example**: Carbon Nanotube (CNT)  
- Purpose: Learn how to set up and analyze SED for **one-dimensional materials**.  

[CNT Example Folder](https://github.com/Tingliangstu/pySED/tree/main/example/CNT)  
<p align="center">
  <img src="https://github.com/Tingliangstu/pySED/blob/main/example/CNT/SED/CNT-SED.svg" alt="CNT SED" width="500">
</p>

### 2D Systems
- **Examples**:  
  - In-plane Graphene
  - [Out-of-plane MoS<sub>2</sub>](https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd)
- Purpose: Learn how to analyze SED for **two-dimensional materials**.  

[Graphene Example Folder](https://github.com/Tingliangstu/pySED/tree/main/example/In_plane_graphene_gpumd)  
<p align="center">
  <img src="https://github.com/Tingliangstu/pySED/blob/main/example/In_plane_graphene_gpumd/SED/compare_LD/Graphene.png" alt="Graphene SED" width="500">
</p>


### 3D Systems
- **Example**: Bulk Silicon  
- Purpose: Learn how to perform SED analysis for **three-dimensional crystalline materials**.  

[Silicon Example Folder](https://github.com/Tingliangstu/pySED/tree/main/example/Silicon_primitive_gpumd)  
<p align="center">
  <img src="https://github.com/Tingliangstu/pySED/blob/main/example/Silicon_primitive_gpumd/SED/compare_LD/Silicon.png" alt="Silicon SED" width="500">
</p>


## Contact info

Ting Liang
liangting.zj@gmail.com;
Wenwu Jiang
wwjiang96@163.com

Or one can raise questions and suggestion at QQ group:
<p align="center">
  <img src="https://github.com/Tingliangstu/pySED/blob/main/docs/source/_static/qq.jpg" alt="Silicon SED" width="400">
</p>