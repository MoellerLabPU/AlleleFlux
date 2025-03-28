Installation
============

Requirements
-----------

AlleleFlux requires Python 3.7 or above and the following dependencies:

* pandas
* numpy
* scipy
* biopython
* tqdm
* pysam
* intervaltree
* statsmodels

For the Snakemake workflow, you'll also need:

* Snakemake 7.0.0 or above

Install from GitHub
------------------

You can install AlleleFlux directly from GitHub:

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/MoellerLabPU/AlleleFlux.git
    cd AlleleFlux
    
    # Install in development mode with PEP 517
    pip install -e . --use-pep517

This will make all the command-line tools available in your environment.

Verify Installation
------------------

To verify that AlleleFlux has been installed correctly, run:

.. code-block:: bash

    alleleflux-profile --help

You should see the help message for the profile command, which confirms that the installation was successful.

Next Steps
---------

Once you've installed AlleleFlux, you can proceed to the :doc:`quickstart` guide to learn how to use it.