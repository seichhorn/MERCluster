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
bootStrapIterations = list(range(int(config['Round1']['bootStrapIterations'])))
# rawDataPath = config['Paths']['rawDataPath']
subclusterTypes = config['Round2']['cellTypes']


rule all:
	input: expand(config['Round2']['Paths']['analysisDir'] + 'cellTypes_{cellType}.txt', cellType = subclusterTypes)

rule full_clustering:
	input:
		rawData

	wildcard_constraints:
		geneSet = '|'.join(config['Round1']['geneSets'])
 
	output:
		fullOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round1']['resolution']), 
		trackingOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round1']['resolution'])

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSet = config['Round1']['geneSets'],
		outputDir = config['Round1']['Paths']['outputDir'],
		clusteringFlavor = config['Round1']['clusteringFlavor'],
		resolution = config['Round1']['resolution']

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor}"


rule bootstrap_clustering:
	input:
		rawData

	output:
		fullOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round1']['resolution']), 
		trackingOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round1']['resolution'])

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSet = config['Round1']['geneSets'],
		outputDir = config['Round1']['Paths']['outputDir'],
		clusteringFlavor = config['Round1']['clusteringFlavor'],
		bootstrapFrac = config['Round1']['bootstrapFrac'],
		resolution = config['Round1']['resolution']

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -fileNameIteration {wildcards.bootstrap} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -bootstrapFrac {params.bootstrapFrac}"

rule select_kValue_R1:
	output:
		config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_8_resolution_1_type_All_geneset_highVar.txt',
		config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_8_resolution_1_type_All_geneset_highVar.txt"

	input:
		expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], geneSet = config['Round1']['geneSets'], cellType = config['Round1']['cellTypes']),
		expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}_bootstrap_{bootstrap}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], geneSet = config['Round1']['geneSets'], cellType = config['Round1']['cellTypes'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round1']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round1']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round1']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"


rule subclustering:
	input:
		experimentData = rawData,
		cellTypes = config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_8_resolution_1_type_All_geneset_highVar.txt',
		cellLabels = config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_8_resolution_1_type_All_geneset_highVar.txt"
	
	wildcard_constraints:
		cellType = '|'.join(subclusterTypes),
		geneSet = '|'.join(config['Round1']['geneSets'])

	output:
		fullOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round2']['resolution']), 
		trackingOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round2']['resolution'])

	params:
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSet = config['Round2']['geneSets'],
		outputDir = config['Round2']['Paths']['outputDir'],
		clusteringFlavor = config['Round2']['clusteringFlavor'],
		resolution = config['Round2']['resolution'],
		restriction = lambda w: config['Round2']['restrictions']["{}".format(w.cellType)]

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.cellType} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -restriction {params.restriction}"


rule subbootstrap_clustering:
	input:
		experimentData = rawData,
		cellTypes = config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_8_resolution_1_type_All_geneset_highVar.txt',
		cellLabels = config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_8_resolution_1_type_All_geneset_highVar.txt"

	wildcard_constraints:
		cellType = '|'.join(subclusterTypes),
		geneSet = '|'.join(config['Round1']['geneSets'])

	output:
		fullOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round2']['resolution']), 
		trackingOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round2']['resolution'])
	
	params:
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSet = config['Round2']['geneSets'],
		outputDir = config['Round2']['Paths']['outputDir'],
		clusteringFlavor = config['Round2']['clusteringFlavor'],
		bootstrapFrac = config['Round2']['bootstrapFrac'],
		resolution = config['Round2']['resolution'],
		restriction = lambda wildcards: config['Round2']['restrictions']['{}'.format(wildcards.cellType)]


	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.cellType} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSet} -fileNameIteration {wildcards.bootstrap} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -restriction {params.restriction} -bootstrapFrac {params.bootstrapFrac}"


rule select_kValue_R2:
	output:
		config['Round2']['Paths']['analysisDir'] + 'cellTypes_{selectedVal}.txt',
		config['Round2']['Paths']['analysisDir'] +"selected_stability_analysis_{selectedVal}.txt"

	input:
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], geneSet = config['Round2']['geneSets'], cellType = config['Round2']['cellTypes']),
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}_bootstrap_{bootstrap}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], geneSet = config['Round2']['geneSets'], cellType = config['Round2']['cellTypes'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round2']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round2']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round2']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"




