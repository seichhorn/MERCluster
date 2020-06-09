import importlib
import networkx

from mercluster.core import analysistask
from mercluster.core import metadataset


class SnakemakeRule(object):

    def __init__(self, analysisTask: analysistask.analysisTask,
                 pythonPath=None):
        self._analysisTask = analysisTask
        self._pythonPath = pythonPath

    @staticmethod
    def _add_quotes(stringIn):
        return '\'%s\'' % stringIn

    @staticmethod
    def _clean_string(stringIn):
        return stringIn.replace('\\', '/')

    def _expand_as_string(self, task, indexCount) -> str:
        fileName = task.analysisName + '_{g}'
        return 'expand(%s, g=list(range(%i)))' % (self._add_quotes(
            task.metaDataSet.get_analysis_path(analysisTask=task,
                                               subDir='tasks',
                                               fileName=fileName,
                                               extension='done')),
                                                  indexCount)

    def _generate_input_names(self, task):
        if task.fragment_count()>1:
            return self._clean_string(self._expand_as_string(
                task, task.fragment_count()))
        else:
            return self._clean_string(self._add_quotes(
                task.metaDataSet.get_analysis_path(analysisTask=task,
                                                   subDir='tasks',
                                                   fileName=task.analysisName,
                                                   extension='done')))

    def _generate_input(self) -> str:
        if len(self._analysisTask.get_dependencies()) > 0:
            if isinstance(self._analysisTask.get_dependencies()[0],
                          analysistask.analysisTask):
                inputTasks = self._analysisTask.get_dependencies()
            else:
                inputTasks = [self._analysisTask.metaDataSet.load_analysis_task(x)
                              for x in self._analysisTask.get_dependencies()]
            inputString = ','.join([self._generate_input_names(x)
                                    for x in inputTasks])
        else:
            inputString = ''
        return self._clean_string(inputString)

    def _generate_output(self) -> str:
        if self._analysisTask.fragment_count()>1:
            return self._clean_string(self._expand_as_string(
                self._analysisTask, self._analysisTask.fragment_count()))
        else:
            return self._clean_string(self._add_quotes(
                self._analysisTask.metaDataSet.get_analysis_path(
                    analysisTask=self._analysisTask,
                    subDir='tasks',
                    fileName=self._analysisTask.analysisName,
                    extension='.done')))

    def _generate_message(self) -> str:
        messageString = \
            ''.join(['Running ', self._analysisTask.analysisName])

        if self._analysisTask.fragment_count()>1:
            messageString += ' {wildcards.i}'

        return self._add_quotes(messageString)

    def _generate_shell(self) -> str:
        if self._pythonPath is None:
            shellString = 'python '
        else:
            shellString = self._clean_string(self._pythonPath) + ' '
        shellString += '-m mercluster '
        shellString += self._clean_string(
            self._analysisTask.metaDataSet.metaDataSetName) + ' '

        shellString += ''.join(
            ['-t ', '\"', self._clean_string(self._analysisTask.analysisName),
             '\"', ' -s \"',
             self._clean_string(self._analysisTask.metaDataSet.analysisHome),
             '\"'])

        if self._analysisTask.fragment_count()>1:
            shellString += ' -i {wildcards.i}'

        return self._add_quotes(shellString)

    def as_string(self) -> str:
        fullString = ('rule %s:\n\tinput: %s\n\toutput: %s\n\tmessage: %s\n\t'
                      + 'shell: %s\n') \
                      % (self._analysisTask.analysisName,
                         self._generate_input(), self._generate_output(),
                         self._generate_message(),  self._generate_shell())

        return fullString

    def full_output(self) -> str:
        return self._generate_input_names(self._analysisTask)


class SnakefileGenerator(object):

    def __init__(self, analysisParameters,
                 metaDataSet: metadataset.metaDataSet,
                 pythonPath: str=None):
        self._analysisParameters = analysisParameters
        self._metaDataSet = metaDataSet
        self._pythonPath = pythonPath

    def _parse_parameters(self):
        analysisTasks = {}
        for tDict in self._analysisParameters['analysis_tasks']:
            analysisModule = importlib.import_module(tDict['module'])
            analysisClass = getattr(analysisModule, tDict['task'])
            analysisParameters = tDict.get('parameters')
            analysisName = tDict.get('analysis_name')
            newTask = analysisClass(
                    self._metaDataSet, analysisParameters, analysisName)
            if newTask.analysisName in analysisTasks:
                raise Exception('Analysis tasks must have unique names. ' +
                                newTask.analysisName + ' is redundant.')
            # TODO This should be more careful to not overwrite an existing
            # analysis task that has already been run.
            newTask.saveTask()
            analysisTasks[newTask.analysisName] = newTask
        return analysisTasks

    def _identify_terminal_tasks(self, analysisTasks):
        taskGraph = networkx.DiGraph()
        for x in analysisTasks.keys():
            taskGraph.add_node(x)

        for x, a in analysisTasks.items():
            for d in a.get_dependencies():
                taskGraph.add_edge(d, x)

        return [k for k, v in taskGraph.out_degree if v == 0
                and not analysisTasks[k].event_status('done')]

    def generate_workflow(self) -> str:
        """Generate a snakemake workflow for the analysis parameters
        of this SnakemakeGenerator and save the workflow into the dataset.

        Returns:
            the path to the generated snakemake workflow
        """
        analysisTasks = self._parse_parameters()
        ruleList = {k: SnakemakeRule(v, self._pythonPath)
                    for k, v in analysisTasks.items()
                    if not v.event_status('done')}

        terminalTasks = self._identify_terminal_tasks(analysisTasks)
        workflowString = 'rule all: \n\tinput: ' + \
            ','.join([ruleList[x].full_output() for x in terminalTasks]) + \
            '\n\n'
        workflowString += '\n'.join([x.as_string() for x in ruleList.values()])

        return self._metaDataSet.save_workflow(workflowString)
