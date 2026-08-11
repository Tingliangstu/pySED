"""Package metadata for pySED."""

from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent
VERSION = {}
exec((ROOT / "pySED" / "version.py").read_text(encoding="utf-8"), VERSION)


setup(
    name="pySED",
    version=VERSION["__version__"],
    description=(
        "Spectral energy density analysis for molecular-dynamics trajectories"
    ),
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Ting Liang",
    author_email="liangting.zj@gmail.com",
    url="https://github.com/Tingliangstu/pySED",
    project_urls={
        "Documentation": "https://pysed.readthedocs.io/",
        "Source": "https://github.com/Tingliangstu/pySED",
        "Issues": "https://github.com/Tingliangstu/pySED/issues",
    },
    license="GPL-3.0-or-later",
    license_files=["LICENSE"],
    packages=find_packages(exclude=("tests", "tests.*")),
    python_requires=">=3.10",
    # CuPy remains optional because its wheel must match the user's CUDA
    # runtime; it is intentionally not installed as a generic dependency.
    install_requires=[
        "numpy",
        "scipy",
        "netCDF4",
        "matplotlib",
        "seaborn",
    ],
    extras_require={
        "test": ["pytest>=7"],
    },
    entry_points={
        "console_scripts": [
            "pySED=pySED.main:main",
            "pysed=pySED.main:main",
            "PYSED=pySED.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
