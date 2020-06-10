from abc import abstractmethod
from mercluster.core import analysistask


class Executor(object):

	def __init__(self):
		super().__init__()

	@abstractmethod
	def run(self, task: analysistask.analysisTask, index: int = None,
			rerunCompleted: bool = False) -> None:
		"""Run an analysis task.

		This method will not run analysis tasks that are already currently
		running and analysis is terminated early due to error or otherwise
		will not be restarted.

		Args:
			task: the analysis task to run.
			index: index of the analysis to run for a parallel analysis task.
			rerunCompleted: flag indicating if previous analysis should be
				run again even if it has previously completed. If overwrite
				is True, analysis will be run on the task regardless of its
				status. If overwrite is False, analysis will only be run on
				the task or fragments of the task that have either not been
				started or have previously completed in error.
		"""
		pass


class LocalExecutor(Executor):

	def __init__(self):
		super().__init__()

	def run(self, task: analysistask.analysisTask, index: int = None,
			rerunCompleted: bool = False) -> None:
		if task.event_status('done', i=index) and not rerunCompleted:
			return

		if rerunCompleted:
			task.reset_analysis(i=index)

		if index is not None:
			task.run(index)
		else:
			task.run()

