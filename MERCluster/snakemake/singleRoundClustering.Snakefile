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


rule all:
	input: config['Round1']['Paths']['analysisDir'] + 'cellTypes.txt'


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
		config['Round1']['Paths']['analysisDir'] + 'cellTypes.txt',
		# config['Round1']['Paths']['analysisDir'] + 'stability_analysis.txt'

	input:
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], geneSet = config['Round1']['geneSets'], cellType = config['Round1']['cellTypes']),
		lambda wildcards: expand(config['Round1']['Paths']['outputDir'] + 'clustering/kValue_{kValue}_resolution_{resolution}_type_{cellType}_geneset_{geneSet}_bootstrap_{bootstrap}.txt', kValue = config['Round1']['kValues'], resolution = config['Round1']['resolution'], geneSet = config['Round1']['geneSets'], cellType = config['Round1']['cellTypes'], bootstrap = bootStrapIterations)

	params:
		inputDir = config['Round1']['Paths']['outputDir'] + 'clustering/',
		outputDir = config['Round1']['Paths']['analysisDir'],
		experimentData = rawData,
		geneIdentityFile = config['Round1']['geneIdentityFile'],

	shell:
		pythonPath+" "+codePath+"MERCluster/example_scripts/defineTypes.py {params.inputDir} {params.outputDir} -experimentData {params.experimentData} -geneIdentityFile {params.geneIdentityFile}"


