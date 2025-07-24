[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)
![pySED Logo](https://github.com/Tingliangstu/pySED/blob/main/docs/source/_static/logo.png)

# pySED
--------------------------
**To implement the [phonon SED (spectral energy desity) method in 2010](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), phonon lifetime can be easily calculated.** 

## Main features
-------------
- Output a very nice SED image 
- Using Lorentz to fitting the SED peak 
- Fit all peaks at once and output phonon lifetimes, or can fit them individually
- Phonon lifetime can only be obtained qualitatively
- Can be used in parallel and get SED results quickly
- Can interface with [GPUMD](https://github.com/brucefan1983/GPUMD?tab=readme-ov-file) and [LAMMPS](https://www.lammps.org/#gsc.tab=0), capable of capturing quantum dynamics.

# Installation instructions
--------------------------

1) From source code (Highly recommended)
```python
python setup.py install --user --prefix=
```
   or one can use:
```python
pip install .
```
For convenience, you may want to copy (or link) the files inside scripts
folder to a location included in $PATH environment variable

2) From pip install (Not currently recommended)

```python
pip install pySED-phonon
```

## Usage
--------------------------
The [example files](https://github.com/Tingliangstu/pySED/tree/main/example) have so many case provided to show how to use **pySED**. 
If one can reproduce the case in the example, one should use the program.

For guidance on using the latest version of [pySED](https://github.com/Tingliangstu/pySED/tree/main), please refer to this (MoS$_2$ example](https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd) for the time being.
I will update the manual as soon as I have time.

**Online manual** [https://pysed.readthedocs.io](https://pysed.readthedocs.io/en/latest/)
	
## Citations

| Reference             | Cite for what?                    |
| --------------------- | --------------------------------- |
| [1]                   | for any work that used `pySED`    |
| [2]                   | fundamental theory on phonon SED |

## References

[1] Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, Jianbin Xu, [PYSED: A tool for extracting kinetic-energy-weighted phonon dispersion and lifetime from molecular dynamics simulations](https://arxiv.org/abs/2505.00353). arXiv preprint arXiv:2505.00353, (2025).
[2] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, A. J. H. McGaughey, [Predicting phonon dispersion relations and lifetimes from the spectral energy density](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), Phys. Rev. B **81**, 081411 (2010).
	
## Publications using pySED

Please visit the [publications page](https://github.com/Tingliangstu/pySED/tree/main/publications).

## Contact info
---------------------------------------------------------
Ting Liang
liangting.zj@gmail.com;
Wenwu Jiang
wwjiang96@163.com

