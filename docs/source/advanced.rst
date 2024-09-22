.. highlight:: rst

Advanced topics
===============

Use non-analytical corrections (NAC) in phonon band structure
-------------------------------------------------------------

Non-analytical corrections (NAC) can be used in phonoLAMMPS during the plot of the
phononn band structure preview. To use them a file named **BORN** containing the Born charges
anf the dielectric tensor should exist in the working directory. This file is formated in phonopy
format (check phonopy documentation for further information). To activate this options **--use_NAC** flag is used ::

    $ phonolammps in.lammps --dim 2 2 2 --use_NAC -p

Notice that this option does not modify the calculated force constants written in **FORCE_CONSTANTS** file.
It is only used in the phonon band structure plot, therefore it makes no effect without **-p** flag.

Atomic position optimization
----------------------------

The atomic positions of the unit cell can be optimized using **--optimization** flag.
This uses LAMMPS minimize to perform a minimization of atomic forces at constant volume.
The forces tolerance for the convergence criterion is defined by **--force_tol** flag. ::

    $ phonolammps in.lammps --dim 2 2 2 --optimization --force_tol 1e-10 -p

This option performs a simple minimization. If a more sophisticated optimization is required
then use LAMMPS directly. This option is best used to refine the atomic positions of an already
optimized unit cell.

Finite temperature force constants
----------------------------------

Finite temperature force constants can be calculated from molecular dynamics using
dynaphopy software (http://abelcarreras.github.io/DynaPhoPy). This software implements
the normal mode projection technique to obtain the renormalized force constants at finite
temperature based on quasiparticle theory. Phonolammps provide a minimum functionality to
make this process automatic using both LAMMPS and dynaphopy python API.
To use this feature dynaphopy must be installed (for further details check dynaphopy documentation).

In command line script temperature is defined by **-t** flag. By default this value is 0 and usual
2n order force constants are calculated. If the temperature is higher than 0 then a molecular
dynamics (MD) simulation is calculated with LAMMPS using a supercell defined by **--dim** flag.
By default the length of the MD if 20 ps with time step of 0.001 ps and a relaxation time of 5 ps,
but these parameters can be tweaked using **--total_time**, **--relaxation_time** and **--relaxation_time** flags.

example for silicon::

    $ phonolammps in.lammps --dim 2 2 2 -pa 0.0 0.5 0.5 0.5 0.0 0.5 0.5 0.5 0.0 -t 300 -p


to have more control over the simulation and the renormalization precedure you will have to use
the two software separately.

