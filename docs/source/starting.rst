Quick start
===========

For **pySED**, MD simulation can be performed by `GPUMD <https://gpumd.org/index.html>`_ or `LAMMPS <https://www.lammps.org/>`_ to generate the required trajectory and velocity information.

.. Note:: 
    After MD simulation, the method of using **pySED** is the same. We provide two examples here, one of which is strongly recommended to repeat before starting your own case.

For GPUMD user
--------------

For GPUMD user, currently, we recommend that users start by learning how to generate SED data for :math:`\text{MoS}_2` using **pySED**. 
This example is well integrated with the latest version of the software and serves as a comprehensive guide.

You can find the :math:`\text{MoS}_2` example here: 
`MoS2 SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_


And the :math:`\text{MoS}_2` example here: 
`Graphene SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_


And the :math:`\text{MoS}_2` example here: 
`Graphene SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>`_

You also can find the :math:`\text{MoS}_2` tutorial here: 
`MoS2 tutorial <https://github.com/Tingliangstu/pySED/blob/main/example/tutorials/MoS2/SED_MoS2.ipynb>`_

The output speed file format is as follows::

   ensemble       nve
   dump_exyz      10     1
   run            500000

This produces a `dump_exyz <https://gpumd.org/gpumd/input_parameters/dump_exyz.html#dump-exyz>`_. file with coordinates and velocity information.


For LAMMPS user
---------------

For LAMMPS users, we recommend starting with the silicon example of **pySED**.

You can find the silicon example here: 
`Silicon SED Example <https://github.com/Tingliangstu/pySED/tree/main/example/Silicon>`_

The output speed file format is as follows::

   dump            vels  all  custom  ${dt_dump}  vels.dat  id  type  vx  vy  vz
   dump_modify     vels  format  line "%d  %d  %0.8g  %0.8g  %0.8g"
   dump_modify     vels  sort  id
   dump            pos   all  custom  ${dt_dump}  pos.dat   id  type  x  y  z
   dump_modify     pos   format  line  "%d  %d  %0.8g  %0.8g  %0.8g"
   dump_modify     pos   sort  id

   run             2097152 

For LAMMPS, the coordinates (`pos.dat`) and velocity (`vels.dat`) files are separate and need to be specified separately in `input_SED.in`.
