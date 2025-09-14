Examples
========

**pySED** provides a rich set of examples to help users understand how to perform Spectral Energy Density (SED) analysis across different types of materials.
These examples are designed to demonstrate best practices and guide users through the workflow using **pySED** in combination with To use **pySED**, molecular dynamics (MD) simulations must first be performed using either  
`GPUMD <https://gpumd.org/index.html>`_ and `LAMMPS <https://www.lammps.org/>`_ .

Recommended Starting Point
--------------------------

We recommend beginning with the following example:

- **Material**: math:`\mathrm{MoS_2}`
- **Purpose**: Demonstrates how to generate SED data using pySED and GPUMD.
- **Link**: `MoS2 SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_

This example is well-integrated with the latest version of pySED and serves as a comprehensive tutorial for new users.

Example Categories
------------------

The examples are organized by material dimensionality.

1D Systems
~~~~~~~~~~

- **Example**: Carbon Nanotube (CNT)
- **Purpose**: Learn how to set up and analyze SED for one-dimensional materials.
- **Link**: `CNT Example Folder <https://github.com/Tingliangstu/pySED/tree/main/example/CNT>`_

2D Systems
~~~~~~~~~~

- **Examples**:
  - In-plane Graphene (`Graphene Example Folder <https://github.com/Tingliangstu/pySED/tree/main/example/In_plane_graphene_gpumd>`_)
  - Out-of-plane math:`\mathrm{MoS_2}` (`MoS2 Example Folder <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_)

- **Purpose**: Learn how to analyze SED for two-dimensional materials.

3D Systems
~~~~~~~~~~

- **Example**: Bulk Silicon
- **Purpose**: Learn how to perform SED analysis for three-dimensional crystalline materials.
- **Link**: `Silicon Example Folder <https://github.com/Tingliangstu/pySED/tree/main/example/Silicon_primitive_gpumd>`_

More Examples
-------------

Additional examples covering various materials and configurations can be found here:
`pySED Example Library <https://github.com/Tingliangstu/pySED/tree/main/example>`_

.. note::

   All examples are compatible with the latest version of pySED.
   Users are encouraged to explore and adapt them to their own research needs.
