Installation
============

Installing **pySED** from source is **highly recommended**, as it is actively
maintained. Use Python 3.10 or newer for pySED 3.0.
The PyPI project named ``pysed`` is a different text-processing utility;
do not use ``pip install pysed`` to install this SED package.

Download the Source Code
---------------------------

You can download the source code from GitHub:

.. code-block:: bash

    git clone https://github.com/Tingliangstu/pySED.git
    cd pySED
    git switch codex/pysed-3.0-architecture

Alternatively, you may download the package directly from the
`pySED 3.0 development branch <https://github.com/Tingliangstu/pySED/tree/codex/pysed-3.0-architecture>`_.
The ``main`` branch currently contains the 2.4.x release.

Install pySED
----------------

Navigate to the downloaded directory and install using one of the following methods:

**Recommended method (via pip):**

.. code-block:: bash

    pip install .

**Alternative method (via setup script):**

.. code-block:: bash

    python setup.py install --user --prefix=

Add Scripts to PATH
-----------------------

For convenience, you may want to copy or link the files inside the `scripts` folder
to a location included in your **$PATH** environment variable.

.. note::

    Depending on your operating system or environment configuration, if you are using conda environment this step may be done automatically.

Verify Installation
----------------------

To confirm that pySED is installed correctly, run:

.. code-block:: bash

    pySED -h

or

.. code-block:: bash

    pysed -h

This will display the help message, including usage instructions and descriptions of input parameters.

.. tip::

    You can always use the `-h` flag to explore available options and understand how to prepare input files for pySED.
