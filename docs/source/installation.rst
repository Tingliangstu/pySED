Installation
============

**pySED** can be installed either from the source code or via the PyPI repository.  
Installing from source is **highly recommended**, as it is actively maintained and provides the latest updates.

Download the Source Code
---------------------------

You can download the source code from GitHub:

.. code-block:: bash

    git clone https://github.com/Tingliangstu/pySED.git

Alternatively, you may download the package directly from the  
`pySED GitHub repository <https://github.com/Tingliangstu/pySED>`_.

Install pySED
----------------

Navigate to the downloaded directory and install using one of the following methods:

**Recommended method (via pip):**

.. code-block:: bash

    cd pySED
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