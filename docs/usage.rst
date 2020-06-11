Usage principles
----------------

Importing MERCluster
~~~~~~~~~~~~~~~~~~~~

Import MERCluster as:

.. code-block:: none

    import mercluster

Constructing an analysis task file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MERCluster uses snakemake_ to coordinate the automated execution of analyses. To
help the user construct and execute snakemake workflows, MERCluster reads in a set
of analysis tasks and parameters in JSON format and converts this to a .Snakefile.

.. _snakemake: https://snakemake.readthedocs.io/en/stable/

MERCluster looks for an entry called ``analysis_tasks``, and retrieves the list
associated with it. Each entry in the list is an analysis task along with any
parameters the user specifies for that task. (an example can be found in
examples/analysis_tasks.json)

.. code-block:: none

    {
            "analysis_tasks": [
        ...
        ]
    }

An entry in the list for the execution of the
mercluster.analysis.cluster.Clustering task could appear as follows:

.. code-block:: none

    {
        "analysis_tasks": [
            ...
            {
                "task": "Clustering",
                "module": "mercluster.analysis.cluster",
                "parameters": {
                    "file_creation_task": "BypassAnalyzedData",
                    "k_values": [10,12,15,20],
                    "resolutions": [1,1.5,2]
                }
            },
            ...
        ]
    }

This would instruct MERCluster to perform Clustering using the parameters contained
within the ``parameters`` entry of this task.

.. note::
    All analysis tasks must have a unique name. If you wanted to run one
    analysis task with different parameters multiple times within the same
    workflow, you can assign an ``analysis_name`` to the task.

    .. code-block:: none

                ...
                {
                    "task": "Clustering",
                    "module": "mercluster.analysis.cluster",
                    "analysis_name": "Clustering_neurons"
                    "parameters": {
                        "file_creation_task": "BypassAnalyzedData",
                        "k_values": [10,12,15,20],
                        "resolutions": [1,1.5,2]
                        "cell_type": "Neurons"
                    }
                },
                {
                    "task": "Clustering",
                    "module": "mercluster.analysis.cluster",
                    "analysis_name": "Clustering_glia"
                    "parameters": {
                        "file_creation_task": "BypassAnalyzedData",
                        "k_values": [10,12,15,20],
                        "resolutions": [1,1.5,2]
                        "cell_type": "Glia"
                    }
                },
                ...

All tasks that the user wants to run should be included, and additional tasks
can be added at a later time if one wants to add additional analyses to a metadataset.
Snakemake automatically determines what needs to be run, and will not re-run tasks
that have already been completed so there is no need to remove completed tasks from the
analysis tasks file.

To create a .Snakefile for a metadataset named metadataset1 that is composed of dataset1 and dataset2,
run the following command:

.. code-block:: none

    python -m mercluster 'metadataset1' --dataset-list dataset1 dataset2 --generate-only -a /path/to/analysistasks.json

Automated execution of analysis tasks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After constructing a .Snakefile for a metadataset and a set of analysis tasks, the
workflow can be executed by running:

.. code-block:: none

    python -m mercluster 'metadataset1'

.. note::
    The construction and execution of a workflow can be performed in one line:

    .. code-block:: none

        python -m mercluster 'metadataset1' --dataset-list dataset1 dataset2 -a /path/to/analysistasks.json

Adding the ``--snakemake-parameters`` flag to these commands allows you to pass
additional parameters to snakemake by providing a path to a json file containing
them. The most typical parameter to pass would be related to executing jobs
on a HPC. Examples of these for execution on an HPC running Slurm can be found in
examples/snakemake_params.json and examples/clusterconfig.json

Execution of a selected task
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A task can be executed outside of the snakemake workflow if desired, to do so
just provide the ``-t`` flag and the name of the task.

.. code-block:: none

    python -m mercluster 'metadataset1' --dataset-list dataset1 dataset2 -t Clustering -a /path/to/analysistasks.json

