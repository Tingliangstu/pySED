# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

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
        "To implement the SED method in 2010, phonon lifetime can be calculated"
    ),
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Liang Ting",
    author_email="liangting.zj@gmail.com",
    url="https://github.com/AI4MD-lab/pySED",
    project_urls={
        "Documentation": "https://pysed.readthedocs.io/",
        "Source": "https://github.com/AI4MD-lab/pySED",
        "Issues": "https://github.com/AI4MD-lab/pySED/issues",
    },
    license="GPL-3.0-or-later",
    license_files=["LICENSE"],
    packages=find_packages(include=("pySED", "pySED.*")),
    include_package_data=False,
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
