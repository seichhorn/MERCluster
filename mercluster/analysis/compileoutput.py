# Most mercluster users will have data in one of two general formats. The first is
# a fully compiled set of results that is ready for use in clustering/downstream
# analysis. In this case, the BypassAnalyzedData allows them to point to the
# prepared data and easily enter the clustering workflow. The second is a set
# of MERlin analysis results that have not been compiled/aggregated/normalized
# beyond the basic processing within MERlin. In this case, the AggregateData task
# will go through each MERlin analysis folder and aggregate the requested task.
# In this case, the CombineOutputs task of MERlin is the one that will
# most frequently be used as the requested task for aggregation.

import os
import errno
from shutil import copy2
from mercluster.core import analysistask


class FileExtensionUnsupported(Exception):
	pass

class MERLinLoadError(Exception):
	pass

class BypassAnalyzedData(analysistask.analysisTask):
	"""
	A metaanalysis task that copies your desired data file to use with
	subsequent methods. Currently designed for csv files. Provide the full
	path to the desired file as the parameter 'source_file'.
	"""

	def __init__(self, metaDataSet, parameters=None, analysisName=None):
		super().__init__(metaDataSet, parameters, analysisName)

		if 'overwrite' not in self.parameters:
			self.parameters['overwrite'] = False

		# If additional loading methods are added then update this list to
		# reflect supported file types
		self.supported_ext = ['.csv', '.h5ad']

		if not os.path.isfile(self.parameters['source_file']):
			logger = self.get_task_logger(self.analysisName)
			error = FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
									self.parameters['source_file'])
			logger.exception(error)
			self.close_task_logger(logger)
			raise error

		self.ext = os.path.splitext(self.parameters['source_file'])[1]
		if self.ext not in self.supported_ext:
			logger = self.get_task_logger(self.analysisName)
			message = 'The source file extension is not currently supported.\
						Currently only {} file types are supported in {}'.\
				format(', '.join(self.supported_ext), self.analysisName)
			logger.error(message)
			self.close_task_logger(logger)
			raise FileExtensionUnsupported(message)

	def get_dependencies(self):
		return []

	def _run_analysis(self, i: int=None) -> None:
		dst = self.metaDataSet.get_analysis_path(analysisTask=self,
												 subDir = 'output',
												 fileName = 'aggregated_data',
												 extension = self.ext)
		copy2(self.parameters['source_file'], dst)

	def return_exported_data(self , **kwargs):
		if self.ext == '.csv':
			return self.metaDataSet.read_csv_to_dataframe('aggregated_data',
														  analysisTask=self,
														  subDir='output',
														  **kwargs)
		elif self.ext == '.h5ad':
			return self.metaDataSet.read_h5ad_to_anndata('aggregated_data',
														 analysisTask=self,
														 subDir='output',
														 **kwargs)

		else:
			print('No method to load {} is currently supported'.format(
				self.ext))


class AggregateMERlinData(analysistask.analysisTask):
	"""
	A metaanalysis task that aggregates data from multiple MERlin datasets.
	Assumes that all exporting tasks in merlin have a method called
	return_exported_data that retrieves their data so that minimal handling
	of filenames/locations is needed to address the different tasks
	that might be requested.
	"""
	def __init__(self, metaDataSet, parameters=None, analysisName=None):
		super().__init__(metaDataSet, parameters, analysisName)

		if 'overwrite' not in self.parameters:
			self.parameters['overwrite'] = False
		self.ext = '.csv'

	def get_dependencies(self):
		return [self.parameters['task_to_aggregate']]

	def return_exported_data(self, **kwargs):
		return self.metaDataSet.read_csv_to_dataframe('aggregated_data',
													  analysisTask=self,
													  subDir='output',
													  **kwargs)

	def _run_analysis(self, i: int=None) -> None:
		allAnalyses = []
		taskToAggregate = self.parameters['task_to_aggregate']
		for dsName in self.metaDataSet.datasetNames:
			ds = self.metaDataSet.load_MERlin_dataset(dsName)
			tempData = ds.load_analysis_task(
				taskToAggregate).return_exported_data()
			allAnalyses.append(tempData)
		if len(allAnalyses) == len(self.metaDataSet.datasetNames):
			combinedAnalysis = pandas.concat(allAnalyses, 0)
			self.metaDataSet.write_csv_from_dataframe(combinedAnalysis,
													  'aggregated_data',
													  analysisTask=self,
													  subDir=output)
		else:
			logger = self.get_task_logger(self.analysisName)
			message = 'Something went wrong in loading the MERlin datasets,\
					   please double check that the requested analyses are\
					   available as your MERlin.ANALYSIS_HOME directory'
			logger.error(message)
			self.close_task_logger(logger)
			raise MERLinLoadError(message)

