import argparse
import os

def parse_args():
	parser = argparse.ArgumentParser(description = 'process options for clustering')

	parser.add_argument('-configFilePath', type = str, help = 'path to use for writing the config file, use .json extension')
	parser.add_argument('-environment', type = str, help = 'name of conda environment for running scripts, give it a full path to the python executable if your environment is not located in ~/.conda/envs/')
	parser.add_argument('-MERClusterLocation', default = '~/MERCluster/', type = str, help = 'path to MERCluster directory' )
	parser.add_argument('-rawDataPath', type = str, help = 'location of raw data file, either .csv or .h5ad')
	parser.add_argument('-outputName', type = str, help = 'base name to use for naming output files, no extension')

	parser.add_argument('-numberOfRounds', default = 2, type = int, help = 'number of rounds of clustering to perform, counting is one-base')
	parser.add_argument('-outputDirs', nargs = '+', type = str, help = 'path for writing output of each round of clustering, supplied as a space separated list of paths ordered with the output for the first round first, second round second, etc')
	parser.add_argument('-analysisDirs', nargs = '+', type = str, help = 'path for writing analysis results for each round of clustering, supplied as a space separated list of paths ordered with the output for the first round first, second round second, etc')

	parser.add_argument('-filterByBatch', nargs = '+', type = str, help = 'True or False with whether to filter datasets based on an observation named batch, enter a value separated by a space for each round of clustering')
	parser.add_argument('-countsCutoffMinMax', nargs = '+', type = float, help = 'Quantiles to use to filter out cells with low or high counts, e.g. 0.01 0.99, if you have multiple rounds enter a value for each round as min max pairs, e.g. 0.1 0.99 0.0 1.0 for two rounds, filtering out the top and bottom 1 percent of cells in the first round and no additional filtering applied in the second round')
	parser.add_argument('-kValues', nargs = '+', type = int, help = 'k values to use for each round of clustering, if you want to use a different set of k values for your different rounds enter a 0 separating each list, e.g. 4 6 8 10 0 10 20 30, would run 4, 6, 8, and 10 for the first round and 10, 20, 30 for the second round')
	parser.add_argument('-bootStrapIterations', nargs = '+', type = int, help = 'number of rounds of bootstrapping to perform, assumes you want the same number for all rounds but you can enter a unique value for each separated by a space if you want')
	parser.add_argument('-geneSets', nargs = '+', type = str, help = 'path to file containing gene sets to use for clustering, if you want to use multiple gene sets in the same round of clustering then enter them as a comma separated list with no spaces between them, separate gene sets to be used in different rounds by a space')
	parser.add_argument('-geneIdentityFile', nargs = '+', type = str, help = 'path to file containing genes to identify cell types, if you supply one file it will be used for all rounds, otherwise supply two paths separated by a space')
	parser.add_argument('-resolution', nargs = '+', type = int, help = 'resolution values to use for each round of clustering, if you want to use a different set of values for your different rounds enter a 0 separating each list, e.g. 1 0 1 2 3 4, would run 1 for the first round and 1, 2, 3, and 4 for the second round')    
	
	parser.add_argument('-cellTypes', nargs = '+', type = str, help = 'Names of cell types to select and recluster in round 2+ of clustering, supply as a comma separated list for cell types to target in a particular round, put a space between comma separated lists to indicate the types for different rounds')    

	args = parser.parse_args()

	return args


def useZeroSep(numList,groupToSelect):
	count = 0
	currentGroup = []
	for element in numList:
		if int(element) == 0:
			count += 1
			if count -1 == groupToSelect:
				break
			else:
				currentGroup = []
		else:
			currentGroup.append(element)
	return currentGroup


def write_config():
	args = parse_args()
	
	if not os.path.exists(os.path.dirname(args.configFilePath)):
		os.makedirs(os.path.dirname(args.configFilePath))

	if len(args.environment.split('/')) == 1:
		envPath = '~/.conda/envs/{}/bin/python'.format(args.environment)
	else:
		envPath = args.environment

	outputDirs = args.outputDirs
	analysisDirs = args.analysisDirs

	countsCutoffMins = [args.countsCutoffMinMax[x] for x in range(len(args.countsCutoffMinMax)) if x%2 == 0] 
	countsCutoffMaxs = [args.countsCutoffMinMax[x] for x in range(len(args.countsCutoffMinMax)) if x%2 == 1] 

	if len(args.geneSets) > 1:
		geneSets = args.geneSets
	else:
		geneSets = args.geneSets*(args.numberOfRounds)

	geneSetsExpanded = []
	for element in geneSets:
		geneSetsExpanded.append(element.split(','))

	if len(args.geneIdentityFile)>1:
		geneIdentityFiles = args.geneIdentityFile
	else:
		geneIdentityFiles = args.geneIdentityFile*(args.numberOfRounds)


	if len(args.cellTypes) == args.numberOfRounds -1:
		cellTypes = args.cellTypes
	elif len(args.cellTypes) > 1 and len(args.cellTypes) < args.numberOfRounds -1:
		print('You have not formatted the cell types list in an acceptable way')
	elif len(args.cellTypes) == 1 and args.numberOfRounds > 2:
		cellTypes = args.cellTypes * (args.numberOfRounds-1)

	cellTypesExpanded = []
	for element in cellTypes:
		cellTypesExpanded.append(element.split(','))




	with open(args.configFilePath,'w') as configOpen:
		configOpen.write('{\n')
		configOpen.write('\t\"Paths\" :\n')
		configOpen.write('\t{\n')
		configOpen.write('\t\t\"pythonDir\" : \"{}\",\n'.format(envPath))
		configOpen.write('\t\t\"codeDir\" : \"{}\",\n'.format(args.MERClusterLocation))
		configOpen.write('\t\t\"rawData\" : \"{}\",\n'.format(args.rawDataPath))
		configOpen.write('\t\t\"fileName\" : \"{}\",\n'.format(args.outputName))
		configOpen.write('\t},\n')

		for clusterRound in range(args.numberOfRounds):
			configOpen.write('\t\"Round{}\" :\n'.format(clusterRound+1))
			configOpen.write('\t{\n')

			configOpen.write('\t\t\"Paths\" :\n')
			configOpen.write('\t\t{\n')
			configOpen.write('\t\t\t\"outputDir\" : \"{}\",\n'.format(outputDirs[clusterRound]))
			configOpen.write('\t\t\t\"analysisDir\" : \"{}\",\n'.format(analysisDirs[clusterRound]))
			configOpen.write('\t\t},\n')

			configOpen.write('\t\t\"Filters\" :\n')
			configOpen.write('\t\t{\n')
			configOpen.write('\t\t\t\"byBatch\" : \"{}\",\n'.format(args.filterByBatch[clusterRound]))
			configOpen.write('\t\t\t\"countsCutoffMin\" : {},\n'.format(countsCutoffMins[clusterRound]))
			configOpen.write('\t\t\t\"countsCutoffMax\" : {},\n'.format(countsCutoffMaxs[clusterRound]))
			configOpen.write('\t\t},\n')
			
			configOpen.write('\t\t\"kValues\" :\n')
			configOpen.write('\t\t[\n')

			if 0 not in args.kValues:
				[configOpen.write('\t\t\t{},\n'.format(x)) for x in args.kValues]
			else:
				[configOpen.write('\t\t\t{},\n'.format(x)) for x in useZeroSep(args.kValues,clusterRound)]
			configOpen.write('\t\t],\n')

			configOpen.write('\t\t\"bootStrapIterations\" :\n')
			configOpen.write('\t\t[\n')

			if len(args.bootStrapIterations) == args.numberOfRounds+1:
				configOpen.write('\t\t\t{},\n'.format(args.bootStrapIterations[clusterRound]))
			else:
				configOpen.write('\t\t\t{},\n'.format(args.bootStrapIterations[0]))

			configOpen.write('\t\t],\n')
			
			configOpen.write('\t\t\"geneSets\" :\n')
			configOpen.write('\t\t[\n')
			[configOpen.write('\t\t\t\"{}\",\n'.format(x)) for x in geneSetsExpanded[clusterRound]]
			configOpen.write('\t\t],\n')

			configOpen.write('\t\t\"geneIdentityFile\" : \"{}\",\n'.format(geneIdentityFiles[clusterRound]))

			configOpen.write('\t\t\"resolution\" :\n')
			configOpen.write('\t\t[\n')

			if 0 not in args.resolution:
				[configOpen.write('\t\t\t{},\n'.format(x)) for x in args.resolution]
			else:
				[configOpen.write('\t\t\t{},\n'.format(x)) for x in useZeroSep(args.resolution,clusterRound)]
			configOpen.write('\t\t],\n')

			if clusterRound > 0:
				configOpen.write('\t\t\"cellTypes\" :\n')
				configOpen.write('\t\t[\n')
				[configOpen.write('\t\t\t\"{}\",\n'.format(x)) for x in cellTypesExpanded[clusterRound-1]]
				configOpen.write('\t\t],\n')

			configOpen.write('\t},\n')


		configOpen.write('}\n')


if __name__ == '__main__':
	write_config()

