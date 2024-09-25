![pySED Logo](https://github.com/Tingliangstu/pySED/blob/main/docs/source/_static/logo.png)
[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)

# pySED
--------------------------
**To implement the phonon SED (spectral energy desity) method in 2010, phonon lifetime can be calculated.** 

## Main features
-------------
- Output a very nice SED image 
- Using Lorentz to fitting the SED peak 
- Fit all peaks at once and output phonon lifetimes, or can fit them individually
- Simple and easy to use
- Phonon lifetime can only be obtained qualitatively
- Can be used in parallel and get SED results quickly

# Installation instructions
--------------------------

1) From source code
```python
# python setup.py install --user --prefix=
```
For convenience, you may want to copy (or link) the files inside scripts
folder to a location included in $PATH environment variable

2) From pip install (Not currently recommended)

```python
# pip install pySED-phonon
```

## Usage
--------------------------
The example file have so many case provided to show how to use **pySED**. 
If one can reproduce the case in the example, one should use the program.

**Online manual** [https://pysed.readthedocs.io](https://pysed.readthedocs.io/en/latest/)

## References
--------------------------
* J. M. Larkin, Ph.D. thesis, Carnegie Mellon University, 2013
* Larkin, et al., *Phys. Rev. B* **81**, 081411(R) (2010)
* J. E. Turney, E. S. Landry, A. J. H. McGaughey, and C. H. Amon, *Phys. Rev. B* **79**, 064301 (2009).
* A. J. H. McGaughey and M. Kaviany, *Phys. Rev. B* **69**, 094303 (2004).
* [https://github.com/tyst3273/phonon-sed](https://github.com/tyst3273/phonon-sed)

## If you use **pySED**, the following citations are highly recommended:

* [Li, J., Ying, P., Liang, T., Du, Y., Zhou, J., & Zhang, J. (2023)](https://doi.org/10.1039/D2CP05673A)
* [Penghua Ying, Ting Liang, Ke Xu, Jin Zhang, Jianbin Xu, Zheng Zhong, and Zheyong Fan (2023)](https://pubs.acs.org/doi/10.1021/acsami.3c07770)

## Contact info
---------------------------------------------------------
Ting Liang
liangting.zj@gmail.com
Wenwu Jiang
wwjiang96@163.com

