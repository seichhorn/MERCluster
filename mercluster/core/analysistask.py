import os
import shutil
import time
import threading
import numpy as np
import logging
from abc import ABC

from typing import List, Optional
from mercluster.utils import logutils

class AnalysisAlreadyStartedException(Exception):
	pass


class AnalysisAlreadyExistsException(Exception):
	pass

class ParameterInconsistencyException(Exception):
	pass

class analysisTask(ABC):
	"""
	Much of the ideas/core code copied from MERlin, doi:10.5281/zenodo.3758540
	but stripped down to retain only required functionality

 	An abstract class for analysis tasks, creates the basic directory structure,
 	sets up logging, and enables running of task.

	"""
	def __init__(self, metaDataSet, parameters=None, analysisName=None):

		if analysisName is None:
			self.analysisName = type(self).__name__
		else:
			self.analysisName = analysisName
		if parameters is None:
			self.parameters = {}
		else:
			self.parameters = parameters

		self.metaDataSet = metaDataSet
		self.analysisBasePath = self.metaDataSet.get_analysis_path(
			analysisTask=self)
		self.taskDir = self.metaDataSet.get_analysis_path(
			analysisTask=self, subDir='tasks')
		self.logPath = self.metaDataSet.get_analysis_path(
			subDir='log', fileName=self.analysisName, extension='.log')

	def get_dependencies(self) -> List:
		return []

	def fragment_count(self) -> int:
		"""Indicates the number of jobs that are associated with this task"""
		return 1

	def saveTask(self, overwrite=False) -> None:
		self.get_task_logger(self.analysisName)

		if overwrite:
			self.remove_all_task_files()

		self.check_parameters(overwrite=overwrite)
		self._create_directory_structure()
		self._save_task_info()

	def _save_task_info(self) -> None:
		"""
		Writes name, class and parameters of the task
		"""
		taskDict = {'parameters': {**self.parameters}}
		taskDict['analysis_module'] = self.__module__
		taskDict['analysis_type'] = type(self).__name__
		self.metaDataSet.write_json_from_dict(taskDict, 'task',
											  analysisTask=self,
											  subDir='tasks')

	def return_task_info(self) -> None:
		"""
		Writes name, class and parameters of the task
		"""
		return self.metaDataSet.read_json_to_dict()('task',
													analysisTask=self,
													subDir='tasks')

	def _create_directory_structure(self):
		"""
		Creates the main analysis directory and all the subdirectories
		"""
		self._create_analysis_subdirectory('tasks')
		self._create_analysis_subdirectory('output')
		self._create_analysis_subdirectory('log')

	def remove_all_task_files(self) -> None:
		self.get_task_logger(self.analysisName)

		shutil.rmtree(self.analysisBasePath, ignore_errors=True)
		logger = self.get_task_logger(self.analysisName)
		logger.info('Deleting all files associated with task')
		self.close_task_logger(logger)

	def check_parameters(self, overwrite=False) -> None:
		parametersPath = self.metaDataSet.get_analysis_path(
			analysisTask=self, fileName='task', extension='.json')
		if overwrite and os.path.exists(parametersPath):
			os.remove(parametersPath)

		if os.path.exists(parametersPath):
			originalParameters = self.metaDataSet.read_json_to_dict(
				'task', analysisTask=self, subDir='tasks')['parameters']
			failures = []
			for k,v in originalParameters.items():
				if k not in self.parameters:
					message = '{} not in current parameters, was present in' \
							  ' parameters originally used to create task'
					failures.append(message.format(k))
				else:
					if self.parameters[k] != v:
						message = 'value inconsistency for {}, original ' \
								  '= {}, current = {}'
						failures.append(message.format(k, v,
													   self.parameters[k]))
			for k,v in self.parameters.items():
				if k not in originalParameters:
					message = '{} not in original parameters, is present ' \
							  'in current parameters'
					failures.append(message.format(k))
				else:
					if originalParameters[k] != v:
						message = 'value inconsistency for {}, original ' \
								  '{}, current = {}'
						failures.append(message.format(k, originalParameters[k],
													   v))
			if len(failures) > 0:
				logger = self.get_task_logger(self.analysisName)
				logger.error('Conflict in current/original parameters of task')
				logger.error('\n'.join(failures))
				self.close_task_logger(logger)
				raise ParameterInconsistencyException('\n'.join(failures))

	def _create_analysis_subdirectory(self, subDir) -> None:
		subDirName = self.metaDataSet.get_analysis_path(analysisTask=self,
														subDir=subDir)
		os.makedirs(subDirName, exist_ok=True)

	def get_task_logger(self, loggerName, i=None,
						overwrite=False) -> logging.Logger:
		if i is None:
			logPath = self.logPath
		else:
			jobName = '{}_{}'.format(self.analysisName, i)
			logPath = self.metaDataSet.get_analysis_path(
				subDir='log', fileName=jobName, extension='.log')
		if overwrite:
			logutils.restart_log(logPath)
		logger = logutils.getLogger(loggerName, logPath)
		return logger

	def close_task_logger(self, logger) -> None:
		logutils.closeLoggerHandlers(logger)

	def run(self, i=None, overwrite=False) -> None:
		jobName = '_'.join([x for x in (self.analysisName, i) if x])
		logger = self.get_task_logger(self.analysisName, i=i)
		try:
			if self.event_status('run', i=i, timeLapsed=30):
				logger.warning('Attempted to start {},'
							   'but it is already running'.format(jobName))
				raise AnalysisAlreadyStartedException(
					'Unable to run {} since it is already running'.format(
						jobName))
			if overwrite:
				self.reset_analysis(i=i)

			if self.event_status('done', i=i):
				logger.warning('Attempted to start {},'
							   'but it is already complete'.format(jobName))
				raise AnalysisAlreadyExistsException(
					'Unable to run {} since it is complete and overwriting was\
					not enabled'.format(jobName))

			logger.info('Starting {}'.format(jobName))
			self._indicate_running(i=i)
			self._run_analysis(i)
			self._record_analysis_event('done', i=i)
			logger.info('Completed {}'.format(jobName))
			self.close_task_logger(logger)
		except Exception as e:
			logger.exception(e)
			self._record_analysis_event('error', i=i)
			self.close_task_logger(logger)
			raise e

	# @abstractmethod
	def _run_analysis(self, i: int=None) -> None:
		"""Method to actually perform actual analysis, implemented by each
		specific subclass"""
		pass

	def _indicate_running(self, i: Optional[int]=None,
						  interval: int=30) -> None:
		"""A loop that regularly signals to the dataset that this analysis
		task is still running successfully.

		Once this function is called, the dataset will be notified every
		time the interval passes that this analysis is still running until
		the analysis completes.
		"""
		if self.event_status('done', i=i) or self.event_status('error',
															   i=i):
			return

		self._record_analysis_event('run', i=i)
		self.runTimer = threading.Timer(interval, self._indicate_running)
		self.runTimer.daemon = True
		self.runTimer.start()

	def _record_analysis_event(self, event: str,
							   i: Optional[int]=None) -> None:
		"""
		Creates a file within the task subdirectory for the current
		analysisTask, the contents is the current time
		Args:
			jobName: name of analysis task
			event: type of event to record (can be running, complete, or error)
		"""
		jobName = '_'.join([x for x in (self.analysisName, i) if x])
		runTaskPath = self.metaDataSet.get_analysis_path(analysisTask=self,
														 subDir='tasks',
														 fileName = jobName,
														 extension = event)
		with open(runTaskPath,'w') as fp:
			fp.write('%s' % time.time())

	def reset_analysis(self, i: Optional[int]=None) -> None:
		"""This clears all job-related files from the task subdirectory,
		it does not remove any other outputs associated with the task
		Args:
			jobName: name of job that tasks should be reset
		"""
		jobName = '_'.join([x for x in (self.analysisName, i) if x])
		allTaskFiles = os.listdir(self.taskDir)
		toRemove = [x for x in allTaskFiles if '{}.'.format(jobName) in x]

		for f in toRemove:
			taskPath = self.metaDataSet.get_analysis_path(analysisTask=self,
														  subDir='tasks',
														  fileName=f)
			os.remove(taskPath)

	def event_status(self, event, i: Optional[int]=None,
					 timeLapsed: int=-1) -> bool:
		"""
		Checks to see if a certain task file exists. If the timeLapsed param
		is passed then it checks to see if the file was created longer than
		that many seconds prior
		Args:
			jobName: Name of task to check, typically the name of the analysis
					 task, if it is a parallel task then with _i appended
					 (e.g., analysisTask, or analysisTask_17)
			event: type of task event to check for, supports running, completed
				   and error
			timeLapsed: time prior to current time that event was allowed to
						occur. Currently just used for running task to try
						and address the case where a task is currently still
						running versus was running a while ago.
		Returns:
			True if file exists and was created at least longer than timeLapsed
			prior, False otherwise
		"""
		if timeLapsed == -1:
			allowedStart = np.NINF
		else:
			allowedStart = time.time() - timeLapsed

		jobName = '_'.join([x for x in (self.analysisName, i) if x])

		if (i is None) and (self.fragment_count() > 1):
			allEvents = []
			for frag in range(self.fragment_count()):
				allEvents.append(self.event_status(event,
												   i=frag,
												   timeLapsed=timeLapsed))
				if False in allEvents:
					return False
				else:
					return True

		else :
			taskFile = self.metaDataSet.get_analysis_path(analysisTask=self,
														  subDir='tasks',
														  fileName=jobName,
														  extension=event)
			if os.path.exists(taskFile):
				with open(taskFile,'r') as fp:
					timeOfEvent = fp.readline()
				if float(timeOfEvent) > allowedStart:
					return True
			else:
				return False
