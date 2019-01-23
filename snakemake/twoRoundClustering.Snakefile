#General snakemake pipeline for scanpy clustering analysis of scseq data
#Author: Stephen Eichhorn
#Email: stephen_eichhorn@fas.harvard.edu
#Version = 1.1

#location for python, primarily for running a python from a conda env
pythonPath = config['Paths']['pythonDir']

#location of python scripts
codePath = config['Paths']['codeDir']

#location of raw data
rawData = config['Paths']['rawData']

# #Parameters for clustering
bootStrapIterations = list(range(int(config['Round1']['bootStrapIterations'][0])))
subclusterTypes = config['Round2']['cellTypes']
# rawDataPath = config['Paths']['rawDataPath']

rule full_clustering:
	input:
		rawData

	output:
		fullOutput = config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution,\d+}.txt', 
		trackingOutput = config['Round1']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{kValue}_resolution_{resolution,\d+}.txt'

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSet = config['Round1']['geneSets'],
		outputDir = config['Round1']['Paths']['outputDir']
	shell:
		pythonPath+" "+codePath+"fullClustering.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -kValue {wildcards.kValue} -resolution {wildcards.resolution}"


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
		pythonPath+" "+codePath+"bootstrapSubset.py {input.experimentData} {params.outputDir} -pathToCellList {input.cellTypes} -cellType {wildcards.types} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -fileNameIteration bootstrap_{wildcards.bootstrap} -cellTypeName type_{wildcards.types} -kValue {wildcards.kValue} -resolution {wildcards.resolution}"

rule select_kValue_R2:
	output:
		config['Round2']['Paths']['analysisDir'] + 'cellTypes.txt',
		config['Round2']['Paths']['analysisDir'] + 'stability_analysis.txt'

	input:
		lambda wildcards: expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{types}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], types = subclusterTypes),
		lambda wildcards: expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{types}_bootstrap_{bootstrap}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], types = subclusterTypes, bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round2']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round2']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round2']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"defineTypes2.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"




