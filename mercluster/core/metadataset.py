import os
import json
import pandas as pd
import importlib
import datetime
from typing import List, Union, Optional, Dict
from merlin.core import dataset
import scanpy as sc
import mercluster
from mercluster.utils import logutils
from mercluster.core import analysistask


listOrPath = Union[List, str]
TaskOrName = Union[analysistask.analysisTask, str]

class FileExtensionUnsupported(Exception):
	pass

class DatasetInconsistencyError(Exception):
	pass

class metaDataSet:
	"""
	Much of the ideas/core code copied from MERlin, doi:10.5281/zenodo.3758540
	but stripped down to retain only required functionality

	 Base class for a metadataset

	"""
	def __init__(self, metaDataSetName, analysisDirectory=None):
		self.metaDataSetName = metaDataSetName

		if analysisDirectory is None:
			self.analysisHome = mercluster.ANALYSIS_HOME
			self.analysisPath = os.path.join(self.analysisHome,
											 self.metaDataSetName)
		else:
			self.analysisHome = analysisDirectory
			self.analysisPath = os.path.join(self.analysisHome,
											 self.metaDataSetName)

		self.logDir = os.path.join(self.analysisPath, 'log')
		self.logPath = os.path.join(self.logDir, self.analysisPath + '.log')

		os.makedirs(self.analysisPath, exist_ok=True)
		os.makedirs(self.logDir, exist_ok=True)

	def load_analysis_task(self,
						   analysisTask: str) -> analysistask.analysisTask:
		try:
			taskInfo = self.read_json_to_dict('task',
											  analysisTask=analysisTask,
											  subDir='tasks')
			currentModule = importlib.import_module(
				taskInfo['analysis_module'])
			currentTask = getattr(currentModule, taskInfo['analysis_type'])
			return currentTask(self, analysisName=analysisTask,
							   parameters=taskInfo['parameters'])
		except Exception as e:
			logger = logutils.getLogger(self.metaDataSetName, self.logPath)
			logger.exception(e)
			logutils.closeLoggerHandlers(logger)
			raise e

	def get_analysis_path(self, analysisTask: Optional[TaskOrName]=None,
						  subDir: Optional[str] = None,
						  fileName: Optional[str]=None,
						  extension: Optional[str]=None) -> str:
		"""
		Constructs the path for an analysis directory or file
		Args:
			analysisTask: The analysis task for the file, if applicable
			subDir: The subdirectory for the file, if applicable
			fileName: The name for the file, if applicable
			extension: The file extension, if applicable
		Returns:
			 path, as a string, based on the supplied values, the starting point
			 is the analysis path of the metadataset.
		"""
		if analysisTask is None:
			analysisName = None
		else:
			if type(analysisTask) == str:
				analysisName = analysisTask
			else:
				analysisName = analysisTask.analysisName

		path = self.analysisPath
		if analysisName is not None:
			path = os.path.join(path, analysisName)

		if subDir is not None:
			path = os.path.join(path, subDir)

		if fileName is not None:
			path = os.path.join(path, fileName)

		if extension is not None:
			if extension[0] != '.':
				extension = '.' + extension
			path = path + extension

		return path

	def write_json_from_dict(self, dataToWrite: Dict, fileName: str,
							 analysisTask: Optional[TaskOrName]=None,
							 subDir: str=None) -> None:
		"""Writes a dictionary to a file in JSON format. Writes file in the base
		directory for the analysis task
		Args:
			dataToWrite: dictionary that should be written to file
			fileName: name of file to write
			analysisTask: Analysis task to associate file with, if any
			subDir: subdirectory to write in, if any

		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.json')
		with open(filePath, 'w') as fp:
			json.dump(dataToWrite, fp, indent=4)

	def read_json_to_dict(self, fileName: str,
						  analysisTask: Optional[TaskOrName]=None,
						  subDir: str=None) -> Dict:
		"""Reads a JSON format file in and returns a dictionary

		Args:
			fileName: name of file to read
			analysisTask: Analysis task associated with file, if any
			subDir: subdirectory to look in, if any
		Returns:
			Dictionary based on json file
		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.json')
		with open(filePath, 'r') as fp:
			return json.load(fp)

	def write_csv_from_dataframe(self, dataToWrite: pd.DataFrame, fileName: str,
								 analysisTask: Optional[TaskOrName]=None,
								 subDir: Optional[str]=None, **kwargs) -> None:
		"""
		writes a csv file form a dataframe.
		Args:
			dataToWrite: pandas dataframe to be written to file
			fileName: Name of file to be written, including the extension
			subDir: subdirectory of current analysis task to write file to
			Additional parameters can be passed to the read_csv method
			from pandas as **kwargs
		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.csv')
		dataToWrite.to_csv(filePath, **kwargs)

	def read_csv_to_dataframe(self, fileName: str,
							  analysisTask: Optional[TaskOrName]=None,
							  subDir: Optional[str]=None,
							  **kwargs) -> pd.DataFrame:
		"""
		Loads a csv file as a dataframe.
		Args:
			fileName: Name of file to load, including the extension
			subDir: subdirectory of current analysis task to look for file
			Additional parameters can be passed to the read_csv method
			from pandas as **kwargs
		Returns:
			a pandas dataframe
		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.csv')

		return pd.read_csv(filePath, **kwargs)

	def write_h5ad_from_anndata(self, dataToWrite: sc.AnnData, fileName: str,
								 analysisTask: Optional[TaskOrName]=None,
								 subDir: Optional[str]=None, **kwargs) -> None:
		"""
		writes a csv file form a dataframe.
		Args:
			dataToWrite: pandas dataframe to be written to file
			fileName: Name of file to be written, including the extension
			subDir: subdirectory of current analysis task to write file to
			Additional parameters can be passed to the read_csv method
			from pandas as **kwargs
		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.h5ad')
		dataToWrite.write(filePath, **kwargs)

	def read_h5ad_to_anndata(self, fileName: str,
							  analysisTask: Optional[TaskOrName]=None,
							  subDir: Optional[str]=None,
							  **kwargs) -> pd.DataFrame:
		"""
		Loads a csv file as a dataframe.
		Args:
			fileName: Name of file to load, including the extension
			subDir: subdirectory of current analysis task to look for file
			Additional parameters can be passed to the read_csv method
			from pandas as **kwargs
		Returns:
			a pandas dataframe
		"""
		filePath = self.get_analysis_path(analysisTask=analysisTask,
										  subDir=subDir,
										  fileName=fileName,
										  extension='.h5ad')

		return sc.read_h5ad(filePath, **kwargs)

	def save_workflow(self, workflowString: str) -> str:
		""" Save a snakemake workflow for analysis of this dataset.

		Args:
			workflowString: a string containing the snakemake workflow
				to save

		Returns: the path to the saved workflow
		"""
		snakemakePath = self.get_snakemake_path()
		os.makedirs(snakemakePath, exist_ok=True)

		workflowPath = os.sep.join(
			[snakemakePath, datetime.datetime.now().strftime('%y%m%d_%H%M%S')]) \
					   + '.Snakefile'
		with open(workflowPath, 'w') as outFile:
			outFile.write(workflowString)

		return workflowPath

	def get_snakemake_path(self) -> str:
		"""Get the directory for storing files related to snakemake.

		Returns: the snakemake path as a string
		"""
		return os.sep.join([self.analysisPath, 'snakemake'])

	def get_latest_snakefile(self) -> str:
		"""
		Get the path to the most recently created snakefile

		Returns:
			the .Snakefile path as a string
		"""
		allSnakeFiles = [x for x in os.listdir(self.get_snakemake_path())
						 if '.Snakefile' in x]
		dates = [x.split('.')[0].split('_') for x in allSnakeFiles]
		mostRecent = sorted(dates, key=lambda x: [x[0], x[1]], reverse=True)[0]
		selected = '_'.join(mostRecent) + '.Snakefile'
		return os.path.join(self.get_snakemake_path(), selected)




class metaMERlinDataSet(metaDataSet):
	"""
	Extends the base class to profile support for loading merlin datasets
	"""
	def __init__(self, metaDataSetName, datasets: Optional[listOrPath]=None,
				 analysisDirectory=None):
		super().__init__(metaDataSetName, analysisDirectory)

		self.datasets = self._load_dataset_info(datasets)
		self.datasetNames = list(self.datasets.keys())

	def _load_dataset_info(self, datasets: Optional[listOrPath]=None) -> Dict:
		"""
		Loads information about the datasets that will be used within this
		metadataset. If there is an existing datasetInfo.json file then it loads
		information from that. Otherwise, it accepts two formats. The first, is
		a list of each dataset's name, this is for convenience to handle the
		simplest use cases. The second is a json format file containing the
		dataset name, and the paths to the data home and analysis home
		directories for that dataset.

		Args:
		 datasets: List of dataset names or a path to a json file with the
		 		   dataset names and analysis/data home paths

		Returns:
			Dictionary of dictionaries keyed on dataset names then on dataHome
			and analysisHome for the corresonding dataset.
		"""
		path = self.get_analysis_path(fileName='datasetInfo',
									  extension='.json')
		pathExists = os.path.exists(path)
		if pathExists:
			with open(path,'r') as fp:
				loadedParams = json.load(fp)

		if type(datasets) == list:
			analysisInfo = dict()
			for dataset in datasets:
				analysisInfo[dataset] = {'dataHome': None,
										 'analysisHome': None}
		elif type(datasets) == str:
			filePath, ext = os.path.splitext(datasets)
			if ext.upper() == '.JSON':
				with open(datasets,'r') as fp:
					importedInfo = json.load(fp)
				analysisInfo = dict()
				for k,v in importedInfo.items():
					analysisInfo[k] = dict()
					if 'dataHome' in v:
						analysisInfo[k]['dataHome'] = v['dataHome']
					else:
						analysisInfo[k]['dataHome'] = None
					if 'analysisHome' in v:
						analysisInfo[k]['analysisHome'] = v['analysisHome']
					else:
						analysisInfo[k]['analysisHome'] = None
			else:
				raise FileExtensionUnsupported(
					"Please provide path to a json file containing the"
					"dataset information")

		if pathExists and (datasets is None):
			return loadedParams
		elif pathExists and (datasets is not None):
			if loadedParams == analysisInfo:
				return loadedParams
			else:
				raise DatasetInconsistencyError(
					'The dataset information you provided conflicts with the'
					'information that was previously saved for this metadataset'
				)
		elif (not pathExists) and (datasets is not None):
			self.write_json_from_dict(analysisInfo, 'datasetInfo')
			return analysisInfo
		else:
			raise FileNotFoundError

	def load_MERlin_dataset(self, dataset: str) -> dataset.MERFISHDataSet:
		"""
		Loads a MERlin MERFISHDataSet
		Args:
			dataset: Name of dataset to load
		returns:
			MERlin MERFISHDataSet
		"""
		return dataset.MERFISHDataSet(
			dataset, dataHome=self.datasets[dataset]['dataHome'],
			analysisHome=self.datasets[dataset]['analysisHome'])
