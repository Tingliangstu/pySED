# pySED
--------------------------
**To implement the phonon SED (spectral energy desity) method in 2010, phonon lifetime can be calculated.** 
**Now, pySED can read [GPUMD](https://gpumd.org/) output**

## Main features
-------------
- Output SED image 
- Using Lorentz to fitting the SED peak 
- Fit all peaks at once and output phonon lifetimes, or can fit them individually
- Simple and easy to use
- Phonon lifetime can only be obtained qualitatively

# Installation instructions
--------------------------

1) From source code
```python
# python setup.py install --user --prefix=
```

For convenience, you may want to copy (or link) the files inside scripts
folder to a location included in $PATH environment variable

## Usage
--------------------------
Here three example are provided to show how to use **pySED**. 
If one can reproduce the case in the example, one should use the program.

If I have the time, I might provide a detailed tutorial manual.

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
* [Zeng Jianhui, Liang Ting, Zhang Jingjing, et al. (2023)](https://doi.org/10.1002/smll.202309338)

## Contact info
---------------------------------------------------------
Ting Liang
liangting.zj@gmail.com

