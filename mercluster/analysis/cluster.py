import os
import numpy as np
import pandas as pd
import contextlib
from typing import List, Tuple

from mercluster.core import analysistask
from mercluster.utils import scanpy_helpers
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=150)  # low dpi (dots per inch) yields small inline figures

class MissingAnnotationError(Exception):
    pass

# TODO it might be helpful to add in a feature that asks if you have input
# a metaMERlinDataSet for this task, and then use the various merlin functions
# to retrieve all the genes in your dataset to make sure your data only has
# genes in it/do that type of preprocessing. Little complicated if your input
# data across various datasets didn't have all the same genes for same reason.
class Clustering(analysistask.analysisTask):
    """
    An analysis task that performs graph-based clustering of cells based on the
    provided data, most typically expression measurements. This method assumes
    that you want to use every entry in your input dataset for clustering.
    """
    def __init__(self, metaDataSet, parameters=None, analysisName=None):
        super().__init__(metaDataSet, parameters, analysisName)

        if 'file_creation_task' not in self.parameters:
            self.parameters['file_creation_task'] = None
        if 'prior_clustering' not in self.parameters:
            self.parameters['prior_clustering'] = None
        if 'cell_type' not in self.parameters:
            self.parameters['cell_type'] = 'All'
        if 'k_values' not in self.parameters:
            self.parameters['k_values'] = [12]
        if 'resolutions' not in self.parameters:
            self.parameters['resolutions'] = [1.0]
        if 'use_PCs' not in self.parameters:
            self.parameters['use_PCs'] = True
        if 'cluster_min_size' not in self.parameters:
            self.parameters['cluster_min_size'] = 10
        if 'clustering_algorithm' not in self.parameters:
            self.parameters['clustering_algorithm'] = 'leiden'

        self.cellType = self.parameters['cell_type']
        self.priorClustering = self.parameters['prior_clustering']

    def fragment_count(self):
        return len(self.parameters['k_values']) *\
               len(self.parameters['resolutions'])

    def get_dependencies(self) -> List:
        if self.parameters['file_creation_task'] is None:
            return []
        else:
            return [self.parameters['file_creation_task']]

    def _load_data(self, overwrite=False) -> sc.AnnData:
        """
        Loads data from the file creation task and returns it as an anndata
        object. If the loaded file is a pandas dataframe it converts it to
        anndata first.

        Returns:
             anndata object, obs are cells, var are genes
        """
        path = os.path.exists(metaDataSet.get_analysis_path(analysisTask=self,
                                                            fileName='aData',
                                                            extension='.h5ad'))
        if overwrite:
            with contextlib.suppress(FileNotFoundError):
                os.remove(path)
        if os.path.exists(path):
            aData = self.metaDataSet.read_h5ad_to_anndata('aData',
                                                          analysisTask=self)
            return aData
        else:
            requestedTask = self.parameters['file_creation_task']
            data = self.metaDataSet.load_analysis_task(
                requestedTask).return_exported_data()
            if type(data) == pd.DataFrame:
                aData = sc.AnnData(data.values)
                aData.obs = data.index.values.tolist()
                aData.var = data.columns.values.tolist()
                self.metaDataSet.write_h5ad_from_anndata(aData,
                                                         analysisTask=self,
                                                         fileName='aData')
                return aData
            elif type(data) == sc.AnnData:
                self.metaDataSet.write_h5ad_from_anndata(aData,
                                                         analysisTask=self,
                                                         fileName='aData')
                return data

    def _expand_k_and_resolution(self) -> List:
        """
        Creates a list of all k value and resolution pairs. Used to use a single
        indexer to get the correct pair of values to run.

        Returns:
            List of kvalue, resolution tuples for all combinations of k and r.
        """
        kValues = self.parameters['k_values']
        resolutions = self.parameters['resolutions']
        allPairs = []
        for k in kValues:
            for r in resolutions:
                allPairs.append([k, r])
        return allPairs

    def _cut_cells_if_requested(self, aData: sc.AnnData,
                                overwrite=False) -> sc.AnnData:
        """
        cuts the anndata object to the cell type of interest. This method first
        looks within the current task base directory for "cellsToUse.csv", loads
        the file assuming no header and the first column is index, and then
        cuts the anndata object to just those cells. If there is no such file,
        or if overwrite is set to True, it uses the task associated with the
        prior_clustering parameter to get the ids for the requested cell type,
        creates the cellsToUse.csv file, then returns the cut anndata object

        This is done so that the user can create a file of the cells they want
        regardless of whether they have done the clustering previously. This
        might come up when filtering out cells based on QC metrics, for example.

        Inputs:
            aData: anndata object to be cut to specific cells
        Returns:
            anndata object cut to requested cells
        """

        kwargs = {'index_col': 0}
        primaryCellPath = self.metaDataSet.get_analysis_path(
            analysisTask=self, fileName='cellsToUse', extension='.csv')
        if overwrite:
            with contextlib.suppress(FileNotFoundError):
                os.remove(primaryCellPath)
        if os.path.exists(primaryCellPath):
            selectedCells = self.metaDataSet.read_csv_to_dataframe(
                'cellsToUse', analysisTask=self, **kwargs)
            return aData[selectedCells.values.tolist(),:]
        else:
            if self.priorClustering is None:
                return aData
            else:
                priorClustering = self.metaDataSet.load_analysis_task(
                    self.priorClustering)
                selectedCells = priorClustering.return_selected_cells(
                    self.cellType)
                self.metaDataSet.write_csv_from_dataframe(selectedCells,
                                                          'cellsToUse',
                                                          analysisTask=self)
                return aData[selectedCells.values.tolist(), :]

    def _select_significant_PCs(self, aData) -> sc.AnnData:
        """
        Computes the PCs of an anndata object and identifies the number that
        explain more variance than is explained by the typical first PC of a
        shuffled instance of the data

        Args:
            aData: anndata object to identify significant PCs
        Returns:
            anndata object with PCs computed. The number of PCs surpassing those
            of the randomized data are stored as a class variable named
            pcsToUse
        """
        maxPCs = int(np.min(aData.X.shape)) - 1
        if maxPCs < 100:
            pcsToCalc = maxPCs
        else:
            pcsToCalc = 100
        sc.tl.pca(aData, svd_solver='arpack', n_comps=pcsToCalc)
        randomPCs = PCA(n_components=1, svd_solver='arpack')
        randomVariance = [randomPCs.fit(
            scanpy_helpers.shuffler(aData.X)).explained_variance_[0]
                          for _ in range(10)]

        # Use only PCs that explain more variance than the random dataframe
        pcsToUse = len(aData.uns['pca']['variance']
                       [aData.uns['pca']['variance'] >
                        np.median(randomVariance)])
        self.pcsToUse = pcsToUse
        print('Using {} PCs'.format(pcsToUse))
        return aData

    def _compute_neighbors(self, aData, kValue) -> sc.AnnData:
        """
        Computes the neighbors for use in clustering. After the initial nearest
        neighbor graph is constructed, the connectivities are revised based on
        the jaccard index of the nodes connected by each edge, transforming
        it into a shared nearest neighbor graph
        Args:
            aData: anndata object to use for calculating shared nearest neighbor
                   graph
            kValue: number of neighbors to use in construction of nearest
                    neighbor graph
        Returns:
            anndata object with 'neighbors' added as unstructred data.
            neighbors contains connectivities, which is the sNN graph, all
            other entries in neighbors are the outputs from the scanpy nearest
            neighbor function
        """
        if self.parameters['use_PCs']:
            sc.pp.neighbors(aData, n_neighbors=int(kValue),
                            n_pcs=self.pcsToUse)
        else:
            sc.pp.neighbors(aData, n_neighbors=int(kValue), n_pcs=0)

        aData.uns['neighbors']['connectivities'] =\
            scanpy_helpers.neighbor_graph(
                scanpy_helpers.jaccard_kernel,
                aData.uns['neighbors']['connectivities'])
        return aData

    def _cluster(self, aData, resolution, clusterMin=10,
                 clusteringAlgorithm='louvain') -> Tuple[pd.DataFrame]:
        """
        Performs the clustering. This function is a little more complicated
        than strictly necessary because it preserves the information about
        the cluster label of each cell during the iterations of the modularity
        optimization. The final result where global modularity has been
        optimized is saved in the task's output subdir, whereas the iteration
        results are saved in output/iterations. It is sometimes useful to expore
        the cluster labels of cells from modularities prior to steady state, as
        they generally reflect coherent groupings that are more granular than
        the final assignments.

        Args:
            aData: anndata object to use for clustering
            resolution: resolution for modularity calculation
            clusterMin: minimum number of cells that must be in a cluster
                        to keep that cluster
            clusteringAlgorithm: choice of algorithm to use for modularity
                                 optimization, currently leiden and louvain are
                                 supported
        Returns:
            a tuple of dataframes, first is a dataframe containig the cluster
            labels from all rounds of modularity optimization, second is just
            the final round of optimization. Index is always cell id
        """
        adjacency = aData.uns['neighbors']['connectivities']
        g = sc.utils.get_igraph_from_adjacency(adjacency, directed=False)

        if clusteringAlgorithm == 'louvain':
            import louvain as clAlgo
            print('using louvain algorithm')
        elif clusteringAlgorithm == 'leiden':
            import leidenalg as clAlgo
            print('using leiden algorithm')

        optimiser = clAlgo.Optimiser()
        tracking = []
        partition = clAlgo.RBConfigurationVertexPartition(
            g, weights='weight', resolution_parameter=resolution)
        partition_agg = partition.aggregate_partition()
        print(partition.summary())

        diff = optimiser.move_nodes(partition_agg)
        while diff > 0.0:
            partition.from_coarse_partition(partition_agg)
            partition_agg = partition_agg.aggregate_partition()
            tracking.append(partition.membership)
            print(partition_agg.summary())
            diff = optimiser.move_nodes(partition_agg)

        df = pd.DataFrame(tracking, columns=aData.obs.index).T

        clusteringOutput = df.iloc[:,[-1]].copy(deep=True)
        colLabel = 'kValue_{}_resolution_{}'.format(self.kValue,
                                                    int(self.resolution))
        clusteringOutput.columns = [colLabel]
        clusteringOutputGrouped = clusteringOutput.groupby(colLabel).size()

        toZero = clusteringOutputGrouped[
            clusteringOutputGrouped < int(min_size)].index.values.tolist()
        mask = clusteringOutput[colLabel].isin(toZero)
        clusteringOutput[colLabel] = clusteringOutput[colLabel].where(~mask,
                                                                      other=-1)
        print('Clustering yields {} clusters with at least {} cells'.format(
            clusteringOutput[colLabel].unique().astype(int).max(), clusterMin))

        return (df, clusteringOutput)

    def _save_clustering_output(self, df : pd.DataFrame, subDir: str,
                                kValue: int, resolution: float,
                                cellType: str) -> None:
        self.metaDataSet.write_csv_from_dataframe(
            df, 'kValue_{}_resolution_{}_type_{}'.format(
                kValue, '_point_'.join(str(resolution).split('.')), cellType),
            analysisTask=self,
            subDir=subDir)

    def return_clustering_result(self, kValue: int, resolution: float,
                                 cellType: str) -> pd.DataFrame:
        """
        Retrieves the cluster labels.

        Args:
            kValue: kvalue of clustering result
            resolution: resolution of clustering result
            cellType: cell type of clustering result

        Returns:
            A pandas dataframe containing the cluster label of each cell.
            A label of -1 indicates a cell that was part of a cluster that was
            below min size.
        """

        data = self.metaDataSet.read_csv_to_dataframe(
            'kValue_{}_resolution_{}_type_{}'.format(
                kValue, '_point_'.join(str(resolution).split('.')), cellType),
            analysisTask=self,
            subDir='output')
        return data

    def return_clustering_iteration_result(self,
                                           kValue,
                                           resolution,
                                           cellType) -> pd.DataFrame:
        """
        Retrieves the full set of cluster labels from each round of modularity
        optimization.

        Args:
            kValue: kvalue of clustering result
            resolution: resolution of clustering result
            cellType: cell types of clustering result
        Returns:
            a pandas dataframe containing the cluster labels from each iteration
            of modularity optimization. Columns are ordered based on rounds, the
            last column is the same as that in the output subdir, but without
            a min cluster size applied.
        """
        data = self.metaDataSet.read_csv_to_dataframe(
            'kValue_{}_resolution_{}_type_{}'.format(kValue, int(resolution),
                                                     cellType),
            analysisTask=self,
            subDir='output/iterations')
        return data

    def _run_analysis(self, i):
        aData = self._load_data()
        kValue, resolution = self._expand_k_and_resolution()[i]
        self.kValue = kValue
        self.resolution = resolution

        aData = self._cut_cells_if_requested(aData)
        sc.pp.scale(aData, max_value=4)
        if self.parameters['use_PCs']:
            aData = self._select_significant_PCs(aData)

        aData = self._compute_neighbors(aData, kValue)

        clusterMin = self.parameters['cluster_min_size']
        clusteringAlgorithm = self.parameters['clustering_algorithm']
        iterations, final = self._cluster(
            aData, resolution, clusterMin=clusterMin,
            clusteringAlgorithm=clusteringAlgorithm)
        self._save_clustering_output(iterations, 'output/iterations',
                                     self.kValue, self.resolution,
                                     self.cellType)
        self._save_clustering_output(final, 'output', self.kValue,
                                     self.resolution, self.cellType)


class BootstrapClustering(Clustering):
    """
    Task to support running clustering on randomly downsampled versions
    of a cell x gene matrix
    """
    def __init__(self, metaDataSet, parameters=None, analysisName=None):
        super().__init__(metaDataSet, parameters, analysisName)

        if 'bootstrap_fraction' not in self.parameters:
            self.parameters['bootstrap_fraction'] = 0.8
        if 'bootstraps' not in self.parameters:
            self.parameters['bootstraps'] = 20
        if 'cluster_task' not in self.parameters:
            self.parameters['cluster_task'] = None

        if self.parameters['cluster_task']:
            clTask = self.metaDataSet.load_analysis_task(
                self.parameters['cluster_task'])
            self.parameters = {**self.parameters, **clTask.parameters}

    def fragment_count(self):
        return len(self.parameters['k_values']) *\
               len(self.parameters['resolutions']) *\
               self.parameters['bootstraps']

    def _expand_k_and_resolution(self):
        kValues = self.parameters['k_values']
        resolutions = self.parameters['resolutions']
        bootstraps = self.parameters['bootstraps']
        allPairs = []
        for b in range(bootstraps):
            for k in kValues:
                for r in resolutions:
                    allPairs.append([k, r, b])
        return allPairs

    def _bootstrap_cells(self, aData: sc.AnnData, frac: float) -> sc.AnnData:
        """
        Takes a random sample from a cell x gene anndata object, taking the
        fraction specified in the task parameters.

        Args:
            aData: Anndata object to be used in downsampling
            frac: Fraction of original cells to keep in downsampling
        Returns:
            a random downsampling of the rows (cells) of an anndata object
        """

        size = int(np.floor(
            aData.shape[0] * frac))
        aData = aData[np.random.choice(aData.obs.index.values,
                                       size=size, replace=False), :].copy()
        return aData

    def _save_clustering_output(self, df : pd.DataFrame, subDir: str,
                                kValue: int, resolution: float,
                                cellType: str, bootstrapNum: int) -> None:
        self.metaDataSet.write_csv_from_dataframe(
            df, 'kValue_{}_resolution_{}_type_{}_bootstrap_{}'.format(
                kValue, '_point_'.join(str(resolution).split('.')), cellType,
                bootstrapNum),
            analysisTask=self,
            subDir=subDir)

    def return_clustering_result(self, kValue: int,
                                 resolution: float, cellType: str,
                                 bootstrapNum: int) -> pd.DataFrame:
        """
        Retrieves the cluster labels.

        Args:
            kValue: kvalue of clustering result
            resolution: resolution of clustering result
            cellType: cell type of clustering result

        Returns:
            A pandas dataframe containing the cluster label of each cell.
            A label of -1 indicates a cell that was part of a cluster that was
            below min size.
        """

        data = self.metaDataSet.read_csv_to_dataframe(
            'kValue_{}_resolution_{}_type_{}_bootstrap_{}'.format(
                kValue, '_point_'.join(str(resolution).split('.')), cellType,
                bootstrapNum),
            analysisTask=self,
            subDir='output')
        return data

    def return_clustering_iteration_result(self, kValue: int,
                                           resolution: float, cellType: str,
                                           bootstrapNum: int) -> pd.DataFrame:
        """
        Retrieves the full set of cluster labels from each round of modularity
        optimization.

        Args:
            kValue: kvalue of clustering result
            resolution: resolution of clustering result
            cellType: cell types of clustering result
        Returns:
            a pandas dataframe containing the cluster labels from each iteration
            of modularity optimization. Columns are ordered based on rounds, the
            last column is the same as that in the output subdir, but without
            a min cluster size applied.
        """
        data = self.metaDataSet.read_csv_to_dataframe(
            'kValue_{}_resolution_{}_type_{}_bootstrap_{}'.format(
                kValue, '_point_'.join(str(resolution).split('.')), cellType,
                bootstrapNum),
            analysisTask=self,
            subDir='output/iterations')
        return data

    def _run_analysis(self, fragmentIndex):
        aData = self._load_data()
        kValue, resolution, bootstrapNum = self._expand_k_and_resolution()[
            fragmentIndex]
        self.kValue = kValue
        self.resolution = resolution
        self.bootstrapNum = bootstrapNum
        aData = self._cut_cells_if_requested(aData)
        aData = self._bootstrap_cells(aData,
                                      self.parameters['bootstrap_fraction'])
        sc.pp.scale(aData, max_value=4)
        if self.parameters['use_PCs']:
            aData = self._select_significant_PCs(aData)

        aData = self._compute_neighbors(aData, kValue)

        clusterMin = self.parameters['cluster_min_size']
        clusteringAlgorithm = self.parameters['clustering_algorithm']
        iterations, final = self._cluster(
            aData, resolution, clusterMin=clusterMin,
            clusteringAlgorithm=clusteringAlgorithm)
        self._save_clustering_output(iterations, 'output/iterations',
                                     self.kValue, self.resolution,
                                     self.cellType, self.bootstrapNum)
        self._save_clustering_output(final, 'output', self.kValue,
                                     self.resolution, self.cellType,
                                     self.bootstrapNum)


class ClusterStabilityAnalysis(analysistask.analysisTask):
    """
    A task that determines the stability of clusters based on
    the proportion of cells originally assigned to a given cluster that
    remain clustered when a random subset of the data is clustered
    """
    def __init__(self, metaDataSet, parameters=None, analysisName=None):
        super().__init__(metaDataSet, parameters, analysisName)

        if 'min_fraction_cells' not in self.parameters:
            self.parameters['min_fraction_cells'] = 0.9
        self.clusterTask = self.metaDataSet.load_analysis_task(
            self.parameters['cluster_task'])
        self.bootstrapTask = self.metaDataSet.load_analysis_task(
            self.parameters['bootstrap_cluster_task'])

        if 'kValues_to_consider' not in self.parameters:
            self.parameters['kValues_to_consider'] = \
                self.clusterTask.parameters['k_values']
        if 'resolutions_to_consider' not in self.parameters:
            self.parameters['resolutions_to_consider'] = \
                self.clusterTask.parameters['resolutions']
        self.parameters['cell_type'] = self.clusterTask.parameters['cell_type']
        self.parameters['bootstrap_iterations'] = \
            self.bootstrapTask.parameters['bootstrap_iterations']

    def get_dependencies(self):
        return [self.parameters['cluster_task'],
                self.parameters['bootstrap_cluster_task']]

    def _gather_data(self, kValue, resolution, cellType, bootstrapIterations):
        clTask = self.metaDataSet.load_analysis_task(
            self.parameters['cluster_task'])
        bootTask = self.metaDataSet.load_analysis_task(
            self.parameters['bootstrap_cluster_task'])

        fullClustering = clTask.return_clustering_result(kValue,
                                                         resolution,
                                                         cellType)
        for result in range(bootstrapIterations):
            bootClustering = bootTask.return_clustering_result(kValue,
                                                               resolution,
                                                               cellType, result)
            if result == 0:
                fullBoot = bootClustering.copy(deep=True)
            else:
                fullBoot = pd.concat([fullBoot, bootClustering], axis=1)

        return fullClustering, fullBoot

    def _determine_stability(self, fullClustering, fullBoot):
        for boot in range(fullBoot.shape[1]):
            tempMerge = fullClustering.merge(fullBoot, left_index=True,
                                             right_index=True)
            tempMerge.columns = ['Full', 'Boot']
            tempMerge = tempMerge[tempMerge['Full'] != -1]
            recovery = tempMerge.groupby(['Full', 'Boot']).size().unstack().\
                max(1).div(tempMerge.groupby('Full').size())
            if boot == 0:
                recoveryDF = pd.DataFrame(recovery)
            else:
                recoveryDF = pd.concat([recoveryDF, pd.DataFrame(recovery)],
                                       axis=1)
        stableClusters = recoveryDF[recoveryDF.median(1) > 0.5]\
            .index.values.tolist()
        colName = fullClustering.columns.values.tolist()[0]

        totalCells = fullClustering.shape[0]
        recoveredCells = len(fullClustering[(fullClustering[colName].isin(
            stableClusters)) & (fullClustering[colName] != -1)].index.
                             values.tolist())

        return (stableClusters, recoveryDF,
                recoveredCells, totalCells)

    def _run_analysis(self) -> None:
        toDF = []
        for kValue in self.parameters['kValues_to_consider']:
            for resolution in self.parameters['resolutions_to_consider']:
                fullClustering, fullBoot = self._gather_data(
                    kValue, resolution, self.parameters['cell_type'],
                    self.parameters['bootstrap_iterations'])

                stableClusters, recoveryDF, recoveredCells, totalCells =\
                    self._determine_stability(fullClustering, fullBoot)
                toDF.append([kValue, resolution, len(stableClusters),
                             len(recoveryDF), recoveredCells, totalCells])
        df = pd.DataFrame(toDF, columns=['kValue', 'resolution',
                                         'stable clusters', 'total clusters',
                                         'stable cells', 'total cells'])
        df['fraction stable clusters'] = df['stable clusters'] /\
                                         df['total clusters']

        df['fraction stable cells'] = df['stable cells'] /\
                                      df['total cells']

        selectedEntry = df[df['fraction stable cells'] >=
                             self.parameters['min_fraction_cells']].sort_values(
            by='stable clusters', ascending=False).iloc[0, :]

        selectedK = selectedEntry.loc['kValue']
        selectedR = selectedEntry.loc['resolution']

        self.metaDataSet.write_csv_from_dataframe(df, 'full_stability_analysis',
                                                  analysisTask=self,
                                                  subDir='output')

        self.metaDataSet.write_json_from_dict(
            {'selected_kValue': selectedK, 'selected_resolution': selectedR},
            'selected_k_and_r', analysisTask=self, subDir='output')

    def retrieve_selected_clustering(self):
        selectedVals = self.metaDataSet.read_json_to_dict(
            'selected_k_and_r', analysisTask=self, subDir='output')

        return self.clusterTask.return_clustering_result(
            selectedVals['selected_kValue'],
            selectedVals['selected_resolution'],
            self.parameters['cell_type'])

    def return_selected_cells(self, cellType: str):
        selectedVals = self.metaDataSet.read_json_to_dict(
            'selected_k_and_r', analysisTask=self, subDir='output')

        selectedClustering = self.retrieve_selected_clustering()
        anno = self.metaDataSet.read_csv_to_dataframe('cellTypeAnnotations',
                                                      analysisTask=self)
        selectedClusters = anno[anno['cell_type'] ==
                                cellType]['cluster'].values.tolist()
        colName = 'kValue_{}_resolution_{}'.format(
            selectedVals['selected_kValue'],
            selectedVals['selected_resolution'])
        return selectedClustering[selectedClustering[colName].isin(
            selectedClusters)].index.values.tolist()






