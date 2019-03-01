import argparse
import os
import json

def parse_args():
	parser = argparse.ArgumentParser(description = 'process options for clustering')

	parser.add_argument('-pathForSnakeFile', type = str, help = 'where to write the snakefile that will run the created pipeline')
	parser.add_argument('-fullClustering', action='store_true', help = 'perform full clustering on data')
	parser.add_argument('-bootStrapClustering', action='store_true', help = 'perform bootstrap clustering on data')
	parser.add_argument('-defineTypes', action='store_true', help = 'run type identification analysis on data')
	parser.add_argument('-rounds', default = 1, type = int, help = 'Number of rounds of clustering to perform')
	parser.add_argument('-tSNE', action='store_true', help = 'perform tSNE calculation and save output')
	parser.add_argument('-config', type = str, help = 'path to a valid config file generated using config_writer.py')


	args = parser.parse_args()

	return args

def MERCluster():
	args = parse_args()

	if args.rounds > 2:
		print('clustering with more than 2 rounds is not currently supported')
		break

	if args.pathForSnakeFile:
		if not os.path.exists(os.path.dirname(args.pathForSnakeFile)):
			os.makedirs(os.path.dirname(args.pathForSnakeFile))

	geneSets = []
	with open(args.config) as configOpen:
    	config = json.load(configOpen)
    	for i in range(1,args.rounds+1):
    		geneSets.append(config['Round{}'.format(i)]['geneSets'])
    		
    # for x in geneSet:
    # 	if x == 'None'

		with open(args.pathForSnakeFile,'w') as snakeFileOpen:
			snakeFileOpen.write('#General snakemake pipeline for scanpy clustering analysis of scseq data\n')
			snakeFileOpen.write('#Author: Stephen Eichhorn\n')
			snakeFileOpen.write('#Email: stephen_eichhorn@fas.harvard.edu\n')
			snakeFileOpen.write('#Version = 1.1\n\n')

			snakeFileOpen.write('#location of python to use, i.e. the python executable in your virtual environment\n')
			snakeFileOpen.write('pythonPath = config[\'Paths\'][\'pythonDir\']\n\n')

			snakeFileOpen.write('#location of python scripts\n')
			snakeFileOpen.write('codePath = config[\'Paths\'][\'codeDir\']\n\n')

			snakeFileOpen.write('#location of raw data\n')
			snakeFileOpen.write('rawData = config[\'Paths\'][\'rawData\']\n\n')


			if args.bootStrapClustering:
				for i in range(args.rounds):
					snakeFileOpen.write('#Parameters for clustering\n')
					snakeFileOpen.write('bootStrapIterations{} = list(range(int(config[\'Round{}\'][\'bootStrapIterations\'])))\n'.format)(i,i)


			if args.rounds:
				if args.rounds > 1:
					for i in range(2,rounds+1):
						snakeFileOpen.write('subclusterTypes{} = config[\'Round{}\'][\'cellTypes\']\n\n'.format)(i,i)

			snakeFileOpen.write('rule all:\n\tinput:\n')
			if args.defineTypes:	
				if args.rounds > 1:
					snakeFileOpen.write('\t\texpand(config[\'Round{}\'][\'Paths\'][\'analysisDir\'] + \'cellTypes_{{types}}.txt\', types=subclusterTypes{})\n'.format(args.rounds,args.rounds))

				elif args.rounds == 1:
					snakeFileOpen.write('\t\tconfig[\'Round{}\'][\'Paths\'][\'analysisDir\'] + \'cellTypes.txt\'\n'.format(args.rounds))

			if args.fullClustering and not args.defineTypes:
				if args.rounds > 1:
					print('You are asking for more than one round of clustering but not providing a way to determine cell types for the later rounds, this is not supported')
				else:
					snakeFileOpen.write('\t\texpand(config[\'Round{}\'][\'Paths\'][\'outputDir\'] + \'clustering/kValue_{{kValue}}_resolution_{{resolution}}.txt\', kValue = config[\'Round{}\'][\'kValues\'], resolution = config[\'Round{}\'][\'resolution\'])\n'.format(args.rounds,args.rounds,args.rounds))

			if args.bootStrapClustering and not args.defineTypes:
				if args.rounds > 1:
					print('You are asking for more than one round of clustering but not providing a way to determine cell types for the later rounds, this is not supported')
				else:
					snakeFileOpen.write('\t\texpand(config[\'Round{}\'][\'Paths\'][\'outputDir\'] + \'clustering/kValue_{{kValue}}_resolution_{{resolution}}_bootstrap_{{bootstrap}}.txt\', kValue = config[\'Round{}\'][\'kValues\'], resolution = config[\'Round{}\'][\'resolution\'], bootstrap = bootStrapIterations)\n'.format(args.rounds,args.rounds,args.rounds))

			if args.tSNE:
				snakeFileOpen.write('\t\tconfig[\'Round1\'][\'Paths\'][\'outputDir\'] + \'processed_data/\' + config[\'Paths\'][\'fileName\'] + \'_tSNE.h5ad\'\n')


			if args.fullClustering:
				for r in range(1,args.rounds+1):
					snakeFileOpen.write('rule cluster{}:\n'.format(r))
					snakeFileOpen.write('\tinput:\n')

					if r == 1:
						snakeFileOpen.write('\t\texperimentData = rawData\n\n')
					else:
						snakeFileOpen.write('\t\texperimentData = rawData,\n')
						snakeFileOpen.write('\t\tcellTypes = config[\'Round{}\'][\'Paths\'][\'analysisDir\'] + \'cellTypes.txt\',\n'.format(r-1))
						snakeFileOpen.write('\t\tcellLabels = config[\'Round{}\'][\'Paths\'][\'analysisDir\'] + \'stability_analysis.txt\'\n\n'.format(r-1))

						snakeFileOpen.write('\twildcard_constraints:\n')
						snakeFileOpen.write('\t\ttypes = \'|\'.join(subclusterTypes{})\n\n'.format(r))

					snakeFileOpen.write('\tparams:\n')
					snakeFileOpen.write('\t\tbyBatch = config[\'Round{}\'][\'Filters\'][\'byBatch\'],\n'.format(r))
					snakeFileOpen.write('\t\tcountsCutoffMin = config[\'Round{}\'][\'Filters\'][\'countsCutoffMin\'],\n'.format(r))
					snakeFileOpen.write('\t\tcountsCutoffMax = config[\'Round{}\'][\'Filters\'][\'countsCutoffMax\'],\n'.format(r))
					snakeFileOpen.write('\t\tgeneSets = config[\'Round{}\'][\'geneSets\'],\n'.format(r))
					snakeFileOpen.write('\t\toutputDir = config[\'Round{}\'][\'Paths\'][\'outputDir\']\n\n'.format(r))

					snakeFileOpen.write('\toutput:\n')
					snakeFileOpen.write('\t\tfullOutput = config[\'Round{}\'][\'Paths\'][\'outputDir\'] + \'clustering/kValue_{{kValue}}_resolution_{{resolution}}_type_{{types}}_geneset_{{genesets}}.txt,\n\n'.format(r))



					snakeFileOpen.write('\tshell:\n')

					snakeFileOpen.write('\t\tpythonPath+\" \"+codePath+\"full_cluster.py {{input.experimentData}} {{params.outputDir}} -pathToCellTypes {{input.cellTypes}} -pathToCellLabels {{input.cellLabels}} -cellType {{wildcards.types}} -byBatch {{params.byBatch}} -countsPercentileCutoffs {{params.countsCutoffMin}} {{params.countsCutoffMax}} -preselectedGenesFile {{params.geneSet}} -cellTypeName type_{{wildcards.types}} -kValue {{wildcards.kValue}} -resolution {{wildcards.resolution}}\n\n')




		


rule bootstrap_clustering:
	input:
		rawData

	output:
		fullOutput = config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_bootstrap_{bootstrap}.txt', 
		trackingOutput = config['Round1']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{kValue}_resolution_{resolution}_bootstrap_{bootstrap}.txt'

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSet = config['Round1']['geneSets'],
		outputDir = config['Round1']['Paths']['outputDir']

	shell:
		pythonPath+" "+codePath+"bootstrapClustering.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -fileNameIteration bootstrap_{wildcards.bootstrap} -kValue {wildcards.kValue} -resolution {wildcards.resolution} "

rule tSNE:
	input:
		rawData

	output:
		out = config['Round1']['Paths']['outputDir'] + 'processed_data/' + config['Paths']['fileName'] + '_tSNE.h5ad' 

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSet = config['Round1']['geneSets'],
		outputDir = config['Round1']['Paths']['outputDir'],
		outName = 'processed_data/' + config['Paths']['fileName'] + '_tSNE.h5ad'

	shell:
		pythonPath+" "+codePath+"generate_tSNE.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -outputFileName {params.outName}"

rule select_kValue_R1:
	output:
		config['Round1']['Paths']['analysisDir'] + 'cellTypes.txt',
		config['Round1']['Paths']['analysisDir'] + 'stability_analysis.txt'

	input:
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution']),
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_bootstrap_{bootstrap}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round1']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round1']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round1']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"


rule subclustering:
	input:
		experimentData = rawData,
		cellTypes = config['Round1']['Paths']['analysisDir'] + 'cellTypes.txt',
		cellLabels = config['Round1']['Paths']['analysisDir'] + 'stability_analysis.txt'
	
	wildcard_constraints:
		types = '|'.join(subclusterTypes)

	output:
		fullOutput = config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{types}.txt', 
		trackingOutput = config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{kValue}_resolution_{resolution}_type_{types}.txt'

	params:
		outputDir = config['Round2']['Paths']['outputDir'],
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSet = config['Round2']['geneSets']

	shell:
		pythonPath+" "+codePath+"clusterSubset.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.types} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -cellTypeName type_{wildcards.types} -kValue {wildcards.kValue} -resolution {wildcards.resolution}"


rule subbootstrap_clustering:
	input:
		experimentData = rawData,
		cellTypes = config['Round1']['Paths']['analysisDir'] + 'cellTypes.txt',
		cellLabels = config['Round1']['Paths']['analysisDir'] + 'stability_analysis.txt'

	wildcard_constraints:
		types = '|'.join(subclusterTypes)

	output:
		fullOutput = config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{types}_bootstrap_{bootstrap}.txt', 
		trackingOutput = config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{kValue}_resolution_{resolution}_type_{types}_bootstrap_{bootstrap}.txt'
	
	params:
		outputDir = config['Round2']['Paths']['outputDir'],
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSet = config['Round2']['geneSets']

	shell:
		pythonPath+" "+codePath+"bootstrapSubset.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.types} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -fileNameIteration bootstrap_{wildcards.bootstrap} -cellTypeName type_{wildcards.types} -kValue {wildcards.kValue} -resolution {wildcards.resolution}"

rule select_kValue_R2:
	output:
		config['Round2']['Paths']['analysisDir'] + 'cellTypes_{types}.txt',
		config['Round2']['Paths']['analysisDir'] + 'stability_analysis_{types}.txt'

	input:
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{{types}}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution']),
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{{types}}_bootstrap_{bootstrap}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round2']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round2']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round2']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile} -cellType {wildcards.types}"






	ex1 = experiment.Experiment(args.dataFile, args.outputLocation)
	ex1.filter(verbose = args.verbose, byBatch = args.byBatch, countsPercentileCutoffs = list(args.countsPercentileCutoffs),genesPercentileCutoffs = list(args.genesPercentileCutoffs), mitoPercentileMax = args.mitoPercentileMax, geneMin = args.geneMin)
	ex1.selectVariableGenes(preselectedGenesFile = args.preselectedGenesFile, dispersionMin = args.dispersionMinMaxThreshold[0], dispersionMax = args.dispersionMinMaxThreshold[1], dispersionThreshold = args.dispersionMinMaxThreshold[2])
	ex1.processData(regressOut=args.regressOut)
	if args.usePCA:
		ex1.selectPCs()
	ex1.computeNeighbors(kValue=args.kValue, usePCA=args.usePCA)
	ex1.cluster(resolution=args.resolution, clusterMin=args.clusterMin, trackIterations=args.trackIterations, fileNameIteration = args.fileNameIteration)

if __name__ == '__main__':
	fullClustering()


