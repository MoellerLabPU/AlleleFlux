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
* matplotlib
* seaborn
* rpy2 (for R integration)

For the Snakemake workflow, you'll also need:

* Snakemake 8.0.0 or above
* snakemake-executor-plugin-cluster-generic

Conda Installation (Recommended)
---------------------------------

The easiest way to install AlleleFlux is using Conda/Mamba with the provided environment file:

.. code-block:: bash

    # Install mamba (if not already installed)
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

    # Clone the repository and install
    git clone https://github.com/MoellerLabPU/AlleleFlux.git
    cd AlleleFlux
    mamba env create -f environment.yml

    # Activate the environment
    mamba activate alleleflux

This will install all dependencies and AlleleFlux in development mode.

Manual Installation from GitHub
-------------------------------

Alternatively, you can install AlleleFlux directly from GitHub:

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/MoellerLabPU/AlleleFlux.git
    cd AlleleFlux
    
    # Install in development mode
    pip install -e .

This will make all the command-line tools available in your environment.

Verify Installation
------------------

To verify that AlleleFlux has been installed correctly, run:

.. code-block:: bash

    # Check if AlleleFlux command-line tools are available
    alleleflux-profile --help
    
Command-line Tools Available
----------------------------

After installation, the following command-line tools are available:

**Analysis scripts:**

* ``alleleflux-profile`` - Profile MAGs using alignment files
* ``alleleflux-allele-freq`` - Analyze allele frequencies  
* ``alleleflux-scores`` - Calculate scores based on allele frequencies
* ``alleleflux-taxa-scores`` - Calculate taxonomic group scores
* ``alleleflux-gene-scores`` - Calculate gene-level scores
* ``alleleflux-outliers`` - Detect outlier genes
* ``alleleflux-cmh-scores`` - Calculate CMH test scores

**Preprocessing scripts:**

* ``alleleflux-metadata`` - Generate MAG metadata
* ``alleleflux-qc`` - Perform quality control
* ``alleleflux-eligibility`` - Generate eligibility tables
* ``alleleflux-preprocess-between-groups`` - Preprocess data between groups
* ``alleleflux-preprocess-within-group`` - Preprocess data within groups

**Statistics scripts:**

* ``alleleflux-lmm`` - Run linear mixed models
* ``alleleflux-single-sample`` - Perform single sample statistical tests
* ``alleleflux-two-sample-paired`` - Perform paired two-sample tests
* ``alleleflux-two-sample-unpaired`` - Perform unpaired two-sample tests
* ``alleleflux-cmh`` - Run Cochran-Mantel-Haenszel tests

**Accessory scripts:**

* ``alleleflux-create-mag-mapping`` - Create MAG mapping files
* ``alleleflux-add-bam-path`` - Add BAM file paths to metadata

Next Steps
---------

Once you've installed AlleleFlux, you can proceed to the :doc:`quickstart` guide to learn how to use it.