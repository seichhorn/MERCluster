import argparse
import cProfile
import json
import sys
import snakemake
import importlib
from typing import TextIO
from typing import Dict
from mercluster.utils import snakewriter


def build_parser():
	parser = argparse.ArgumentParser(description='Analyze datasets',
									 argument_default=None)

	parser.add_argument('--profile', action='store_true',
						help='enable profiling')
	parser.add_argument('--generate-only', action='store_true',
						help='only generate the directory structure and ' +
							 'do not run any analysis.')
	parser.add_argument('metadataset',
						help='name of metadataset')
	parser.add_argument('--metadataset-class',
						default="metaMERlinDataSet",
						help='the type analysis performed on the datasets')

	group = parser.add_mutually_exclusive_group(required=False)
	group.add_argument('--dataset-list')
	group.add_argument('--dataset-path')
	parser.add_argument('--dataset-list', nargs='*',
						help='list of dataset names that comprise the '
							 + 'metadataset when dataset class is '
							 + 'mercluster.core.metadataset.metaMERlinDataSet,'
							 + 'use either --dataset-list or --dataset-path')
	parser.add_argument('--dataset-path',
						help='path to json file specifying dataset, use either '
							 + '--dataset-list or --dataset-path')

	parser.add_argument('-a', '--analysis-parameters',
						help='path to the analysis parameters file to use')
	parser.add_argument(
		'-t', '--analysis-task',
		help='the name of the analysis task to execute. If no '
			 + 'analysis task is provided, all tasks are executed.')
	parser.add_argument(
		'-i', '--fragment-index', type=int,
		help='the index of the fragment of the analysis task to execute')
	parser.add_argument('-s', '--analysis-home',
						help='the analysis home directory')
	parser.add_argument('--snakefile-path',
						help='path to a snakefile, only necessary if running'
							 'with a snakefile other than the one most recently'
							 'created for this metadataset')
	parser.add_argument('-k', '--snakemake-parameters',
						help='the name of the snakemake parameters file')

	return parser

def mercluster():
	print('MERCluster')
	parser = build_parser()
	args, argv = parser.parse_known_args()

	if args.profile:
		profiler = cProfile.Profile()
		profiler.enable()

	if args.dataset_path is not None:
		datasetNames = args.dataset_path
	elif args.dataset_list is not None:
		datasetNames = args.dataset_list
	else:
		datasetNames = None

	metaDataSetArgs = {'datasets': datasetNames,
					   'analysisDirectory': args.analysis_home}

	metaDataSetModuleName = 'mercluster.core.metadataset'
	metaDataSetClassName = args.metadataset_class

	metaDataSetModule = importlib.import_module(metaDataSetModuleName)
	metaDatasetClass = getattr(metaDataSetModule, metaDataSetClassName)
	metaDataSet = metaDatasetClass(args.metadataset, **metaDataSetArgs)

	snakefilePath = None
	if args.analysis_parameters:
		# This is run in all cases that analysis parameters are provided
		# so that new analysis tasks are generated to match the new parameters
		with open(args.analysis_parameters, 'r') as f:
			snakefilePath = generate_analysis_tasks_and_snakefile(metaDataSet,
																  f)

	if not args.generate_only:
		if args.analysis_task:
			print('Running %s' % args.analysis_task)
			e.run(metaDataSet.load_analysis_task(args.analysis_task),
				  index=args.fragment_index)
		else:
			if args.snakemake_parameters:
				with open(args.snakemake_parameters) as f:
					snakemakeParameters = json.load(f)
			else:
				snakemakeParameters = {}
			if args.snakefile_path:
				snakefilePath = args.snakefile_path
			elif snakefilePath is None:
				snakefilePath = metaDataSet.get_latest_snakefile()
			run_with_snakemake(metaDataSet, snakefilePath, snakemakeParameters)

def generate_analysis_tasks_and_snakefile(dataSet,
										  parametersFile: TextIO) -> str:
	print('Generating analysis tasks from %s' % parametersFile.name)
	analysisParameters = json.load(parametersFile)
	snakeGenerator = snakewriter.SnakefileGenerator(
		analysisParameters, dataSet, sys.executable)
	snakefilePath = snakeGenerator.generate_workflow()
	print('Snakefile generated at %s' % snakefilePath)
	return snakefilePath

def run_with_snakemake(
		dataSet, snakefilePath: str,
		snakemakeParameters: Dict = {},):
	print('Running MERlin pipeline through snakemake')
	snakemake.snakemake(snakefilePath, workdir=dataSet.get_snakemake_path(),
						lock=False,	**snakemakeParameters)
