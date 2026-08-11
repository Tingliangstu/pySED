Installation
============

Use an isolated Python 3.10-or-newer environment.

From PyPI
---------

.. code-block:: bash

   python -m pip install pySED

From source
-----------

.. code-block:: bash

   git clone https://github.com/Tingliangstu/pySED.git
   cd pySED
   python -m pip install .

Verify the command:

.. code-block:: bash

   pySED -h

Optional CuPy backend
---------------------

Install the wheel matching the installed CUDA Toolkit, for example:

.. code-block:: bash

   python -m pip install cupy-cuda12x
   # or, for a CUDA 13 installation:
   python -m pip install cupy-cuda13x

Then select ``backend = cupy``. A normal CPU installation does not require
CUDA or CuPy.
