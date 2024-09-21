try:
    from setuptools import setup, Extension, find_packages
    use_setuptools = True
    print("setuptools is used.")
except ImportError:
    from distutils.core import setup, Extension, find_packages
    use_setuptools = False
    print("distutils is used.")


def get_version_number():
    for l in open('pySED/version.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__

setup(name='pySED',
      version=get_version_number(),
      description='SED method in 2010',
      long_description='To implement the SED method in 2010, phonon lifetime can be calculated',
      author='Liang Ting',
      url='https://github.com/Tingliangstu/pySED',
      author_email='liangting.zj@gmail.com',
      packages=find_packages(),
      scripts=['scripts/pySED'],
      python_requires=">=3.7",
      install_requires=[
            "numpy>=1.15.0",
            "matplotlib>=3.5.2",
            "seaborn",
            "h5py",
            "scipy"],
      license='MIT License'
      )