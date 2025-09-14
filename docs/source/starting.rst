Quick Start
===========

To use **pySED**, molecular dynamics (MD) simulations must first be performed using either  
`GPUMD <https://gpumd.org/index.html>`_ or `LAMMPS <https://www.lammps.org/>`_ to generate the required trajectory and velocity data.

.. note::

   After MD simulation, the usage of **pySED** is the same regardless of the MD engine.  
   We provide two examples below, and strongly recommend repeating one before applying pySED to your own system.

GPUMD Users
-----------

For GPUMD users, we recommend starting with the :math:`\text{MoS}_2` example.  
This case is well-integrated with the latest version of pySED and serves as a comprehensive guide.

- Example folder:  
  `:math:`\text{MoS}_2` SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_

- Jupyter tutorial:  
  `:math:`\text{MoS}_2` Tutorial Notebook <https://github.com/Tingliangstu/pySED/blob/main/example/tutorials/MoS2/SED_MoS2.ipynb>`_

**Required GPUMD input format:**

.. code-block:: bash

   ensemble       nve
   dump_exyz      10     1
   run            500000

This produces a `dump_exyz <https://gpumd.org/gpumd/input_parameters/dump_exyz.html#dump-exyz>`_ file containing atomic coordinates and velocities.

LAMMPS Users
------------

For LAMMPS users, we recommend starting with the silicon example.

- Example folder:  
  `Silicon SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/For_old_version_example/Silicon>`_

**Required LAMMPS input format:**

.. code-block:: bash

   dump            vels  all  custom  ${dt_dump}  vels.dat  id  type  vx  vy  vz
   dump_modify     vels  format  line "%d  %d  %0.8g  %0.8g  %0.8g"
   dump_modify     vels  sort  id
   dump            pos   all  custom  ${dt_dump}  pos.dat   id  type  x  y  z
   dump_modify     pos   format  line  "%d  %d  %0.8g  %0.8g  %0.8g"
   dump_modify     pos   sort  id

   run             2097152

In LAMMPS, the coordinate (`pos.dat`) and velocity (`vels.dat`) files are generated separately.  
These files must be specified in the `input_SED.in` file when running pySED.

.. tip::

   Always verify that your trajectory files are correctly formatted before running pySED.  
   You can use the `pySED -h` command to check available options and input requirements.