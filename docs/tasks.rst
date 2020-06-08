Analysis tasks
==============

.. _compileoutput-label:

compileoutput
-------------

BypassAnalyzedData
^^^^^^^^^^^^^^^^^^

| Description:
| Creates a copy of a user-selected .csv file. This enables easy integration
 with any of the downstream analysis tasks

| Parameters:
| * source\_file -- Path to file that should be copied
| * overwrite -- Boolean to indicate whether to overwrite any files present in this analysis task

AggregateMERlinData
^^^^^^^^^^^^^^^^^^^
| Description:
| Loads the output of a selected analysis task from MERlin for each dataset included
 in the metaDataSet and combines them into a single, vertically concatenated file.

| Parameters:
| * task\_to\_aggregate -- Name of the MERlin task for which the analysis result should be concatentated
| * overwrite -- Boolean to indicate whether to overwrite any files present in this analysis task

cluster
-------------------

Clustering
^^^^^^^^^^
| Description:
| performs graph-based clustering of cells based on the provided data, most
 typically expression measurements. This method assumes that you want to use
 every entry in your input dataset for clustering

| Parameters:
| * file\_creation\_task -- Name of the :ref:`compileoutput` task that contains the expression data
| * prior\_clustering -- Name of the :ref:`Clustering` task, if any, that preceded this round of clustering
 (this is typically used when reclustering a subset of cells that were partitioned out in a prior clustering round.
| * cell\_type -- The name of the type of cells that are to be clustered (used in conjunction with prior\_clustering, see examples/cellTypeAnnotations.csv).
| * k\_values -- A list of the k values to use when constructing the k-nearest neighbor graph
| * resolutions -- A list of the resolutions to use for clustering
| * use\_PCs -- Boolean indicating whether to reduce dimensionality with PCA. Only PCs explaining more variance than the 1st PC of a randomized version of the dataset are kept
| * cluster\_min\_size -- Minimum number of cells that must be in a cluster
| * clustering\_algorithm -- Modularity optimization algorithm to use, either leiden or louvain

BootstrapClustering
^^^^^^^^^^^^^^^^^^^
| Description:
| Performs clustering on randomly downsampled (in terms of rows) instance of the input data matrix.

| Parameters:
| * bootstrap\_frac -- Fraction of rows to retain
| * bootstraps -- Number of different downsamplings to analyze for a given k value and resolution pairing
| * cluster\_task -- Name of the :ref:`Clustering` task that this bootstrap analysis is associated with (this
 ensures that the same set of clustering parameters are used between that clustering and the bootstrap clusterings).

ClusterStabilityAnalysis
^^^^^^^^^^^^^^^^^^^^^^^^
| Description:
| Calculates the stability of clusters based on jaccard similarity of the most
 pair of clusters from a full clustering result and a bootstrap clustering result when
 analyzed with the same k value and resolution, for all k value and resolution pairs analyzed.
 The clustering parameters that yielded the clustering result that contains at least 90% of
 the cells (or whatever is set as min\_fraction\_cells) in stable clusters and also yielded
 the greatest number of clusters is selected, and the full stability metrics are output to a table.
 If a cellTypeAnnotations.csv file is placed in the base directory of this analysis task it can be used
 to retrieve a subset of cells based on the labels on that annotation and the selected clustering result.

| Parameters:
| * min\_fraction\_cells: Minimum number of cells that must be in stable clusters for a result to be considered
| * cluster\_task: Name of the :ref:`Clustering` task to get results from
| * bootstrap_cluster_task: name of the :ref:`BootstrapClustering` task to get results from
| * kValues\_to\_consider: list of k values to use in stability analysis if restricting to a subset of those analysed in the clustering
| * resolutions\_to\_consider: list of resolutions to use in stability analysis if restricting to a subset of those analysed in the clustering
