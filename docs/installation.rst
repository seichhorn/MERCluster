Installation
------------

If you are working on a shared system it is best to install MERCluster into a
new virtual environment. MERCluster requires python 3.6 or above.

To create a virtual environment with conda (if you are on a system without
conda, consider installing Miniconda_ first):

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html

.. code-block:: none

    conda create -n mercluster python=3.6

Activate the environment:

.. code-block:: none

    conda activate mercluster

Install MERCluster within the environment you just created and activated:

.. code-block:: none

    git clone https://github.com/seichhorn/MERCluster.git
    pip install -e MERCluster

If you plan to use MERCluster to interface directly with MERlin analyses,
install MERlin as follows:

.. code-block:: none

    conda install rtree pytables
    git clone https://github.com/ZhuangLab/MERlin.git
    pip install -e MERlin

If you wish to use this environment in a jupyter notebook, add it as a kernel:

.. code-block:: none

    conda install -c anaconda ipykernel
    python -m ipykernel install --user --name=mercluster

You can now deactivate the environment:

.. code-block:: none

    conda deactivate

.. note::
    For versions of conda prior to 4.6, ``conda activate`` and ``conda deactivate`` will
    need to be replaced with platform-specific commands:

        * Windows: ``activate`` or ``deactivate``
        * Linux and macOS: ``source activate`` or ``source deactivate``





