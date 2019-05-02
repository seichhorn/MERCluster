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
geneSetNamesList1 = [k for k,v in list(config['Round1']['geneSets'].items())]
geneSetNamesList2 = [k for k,v in list(config['Round2']['geneSets'].items())]


rule all:
	input: expand(config['Round2']['Paths']['analysisDir'] + 'cellTypes_{cellType}_analysis_geneset_{geneSet}.txt', cellType = subclusterTypes, geneSet = geneSetNamesList2)


rule full_clustering:
	input:
		rawData

	wildcard_constraints:
		geneSet = '|'.join(geneSetNamesList1)
 
	output:
		fullOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round1']['resolution']), 
		trackingOutput = expand(config['Round1']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet}}.txt', resolution = config['Round1']['resolution'])

	params:
		byBatch = config['Round1']['Filters']['byBatch'],
		countsCutoffMin = config['Round1']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round1']['Filters']['countsCutoffMax'],
		geneSetPath = lambda w: config['Round1']['geneSets']["{}".format(w.geneSet)],
		outputDir = config['Round1']['Paths']['outputDir'],
		clusteringFlavor = config['Round1']['clusteringFlavor'],
		merfish = config['Round1']['merfish'],
		resolution = config['Round1']['resolution']


	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSetPath} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -merfish {params.merfish}"


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
		geneSetPath = lambda w: config['Round1']['geneSets']["{}".format(w.geneSet)],
		outputDir = config['Round1']['Paths']['outputDir'],
		clusteringFlavor = config['Round1']['clusteringFlavor'],
		bootstrapFrac = config['Round1']['bootstrapFrac'],
		merfish = config['Round1']['merfish'],
		resolution = config['Round1']['resolution']

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input} {params.outputDir} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSetPath} -fileNameIteration {wildcards.bootstrap} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -bootstrapFrac {params.bootstrapFrac} -merfish {params.merfish}"

rule select_kValue_R1:
	output:
		expand(config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_10_resolution_1_type_All_geneset_{geneSet}.txt', geneSet = geneSetNamesList1),
		expand(config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_10_resolution_1_type_All_geneset_{geneSet}.txt", geneSet = geneSetNamesList1)

	input:
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], cellType = config['Round1']['cellTypes'], geneSet = geneSetNamesList1),
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}_bootstrap_{bootstrap}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], cellType = config['Round1']['cellTypes'], geneSet = geneSetNamesList1, bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round1']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round1']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round1']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"


rule subclustering:
	wildcard_constraints:
		cellType = '|'.join(subclusterTypes),
		geneSet1 = '|'.join(geneSetNamesList1),
		geneSet2 = '|'.join(geneSetNamesList2)

	input:
		experimentData = rawData,
		cellTypes = expand(config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_10_resolution_1_type_All_geneset_{geneSet1}.txt', geneSet1 = geneSetNamesList1),
		cellLabels = expand(config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_10_resolution_1_type_All_geneset_{geneSet1}.txt", geneSet1 = geneSetNamesList1)


	output:
		fullOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet2}}.txt', resolution = config['Round2']['resolution']), 
		trackingOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet2}}.txt', resolution = config['Round2']['resolution'])

	params:
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSetPath = lambda w: config['Round2']['geneSets']["{}".format(w.geneSet2)],
		outputDir = config['Round2']['Paths']['outputDir'],
		clusteringFlavor = config['Round2']['clusteringFlavor'],
		merfish = config['Round1']['merfish'],
		resolution = config['Round2']['resolution'],
		restriction = lambda w: config['Round2']['restrictions']["{}".format(w.cellType)]

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.cellType} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSetPath} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -restriction {params.restriction} -merfish {params.merfish}"


rule subbootstrap_clustering:
	wildcard_constraints:
		cellType = '|'.join(subclusterTypes),
		geneSet1 = '|'.join(geneSetNamesList1),
		geneSet2 = '|'.join(geneSetNamesList2)

	input:
		experimentData = rawData,
		cellTypes = expand(config['Round1']['Paths']['analysisDir'] + 'cellTypes_kValue_10_resolution_1_type_All_geneset_{geneSet1}.txt', geneSet1 = geneSetNamesList1),
		cellLabels = expand(config['Round1']['Paths']['analysisDir'] +"selected_stability_analysis_kValue_10_resolution_1_type_All_geneset_{geneSet1}.txt", geneSet1 = geneSetNamesList1)

	output:
		fullOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet2}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round2']['resolution']), 
		trackingOutput = expand(config['Round2']['Paths']['outputDir'] + 'clustering/iterationTracking/kValue_{{kValue}}_resolution_{resolution}_type_{{cellType}}_geneset_{{geneSet2}}_bootstrap_{{bootstrap}}.txt', resolution = config['Round2']['resolution'])
	
	params:
		byBatch = config['Round2']['Filters']['byBatch'],
		countsCutoffMin = config['Round2']['Filters']['countsCutoffMin'],
		countsCutoffMax = config['Round2']['Filters']['countsCutoffMax'],
		geneSetPath = lambda w: config['Round2']['geneSets']["{}".format(w.geneSet2)],
		outputDir = config['Round2']['Paths']['outputDir'],
		clusteringFlavor = config['Round2']['clusteringFlavor'],
		bootstrapFrac = config['Round2']['bootstrapFrac'],
		merfish = config['Round1']['merfish'],
		resolution = config['Round2']['resolution'],
		restriction = lambda wildcards: config['Round2']['restrictions']['{}'.format(wildcards.cellType)]


	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/cluster.py {input.experimentData} {params.outputDir} -pathToCellTypes {input.cellTypes} -pathToCellLabels {input.cellLabels} -cellType {wildcards.cellType} -byBatch {params.byBatch} -countsPercentileCutoffs {params.countsCutoffMin} {params.countsCutoffMax} -preselectedGenesFile {params.geneSetPath} -fileNameIteration {wildcards.bootstrap} -kValue {wildcards.kValue} -resolution {params.resolution} -clusteringAlgorithm {params.clusteringFlavor} -restriction {params.restriction} -bootstrapFrac {params.bootstrapFrac} -merfish {params.merfish}"


rule select_kValue_R2:
	output:
		config['Round2']['Paths']['analysisDir'] + 'cellTypes_{cellType}_analysis_geneset_{geneSet}.txt',

	input:
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], geneSet = geneSetNamesList2, cellType = config['Round2']['cellTypes']),
		expand(config['Round2']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}_bootstrap_{bootstrap}.txt', kValue = config['Round2']['kValues'], resolution = config['Round2']['resolution'], geneSet = geneSetNamesList2, cellType = config['Round2']['cellTypes'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round2']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round2']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round2']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"




