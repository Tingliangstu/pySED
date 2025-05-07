try:
    from setuptools import setup, find_packages
    use_setuptools = True
    print("setuptools is used.")
except ImportError:
    from distutils.core import setup, find_packages
    use_setuptools = False
    print("distutils is used.")

def get_version_number():
    with open('pySED/version.py', 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                exec(line, globals())
                return __version__
    raise ValueError("Version number not found.")

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pySED',
    version=get_version_number(),
    description='To implement the SED method in 2010, phonon lifetime can be calculated',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
        "scipy",
        "dynasor>=2.0",
    ],
    license='MIT License',
)
