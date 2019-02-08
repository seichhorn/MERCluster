import argparse
from MERCluster.classes import experiment

def parse_args():
	parser = argparse.ArgumentParser(description = 'process options for clustering')

	parser.add_argument('dataFile', type = str, help = 'path to counts matrix with rows as cells and columns as genes saved as .h5ad')
	parser.add_argument('outputLocation', type = str, help = 'directory to write results')
	parser.add_argument('-byBatch', default = False, type = bool, help = 'if the batch variable is present in the AnnData object you can set this to true to apply percentile based filters to each dataset independently' )
	parser.add_argument('-countsPercentileCutoffs', nargs = 2, default = [0.0, 1.0], type = float, metavar = ('min counts percentile','max counts percentile'), help = 'cutoff values to remove cells below/above a certain minimum or maximum number of counts, represented as a percentile, enter minimum first then maximum')
	parser.add_argument('-genesPercentileCutoffs', nargs = 2, default = [0.0, 1.0], type = float, metavar = ('min genes percentile','max genes percentile'), help = 'cutoff values to remove cells below/above a certain minimum or maximum number of genes, represented as a percentile, enter minimum first then maximum')
	parser.add_argument('-mitoPercentileMax', default = 0.02, type = float, help = 'cutoff value to remove cells above a maximum percent of mitochondrial reads, set to None if you want to skip this')
	parser.add_argument('-geneMin', default = 10, type = int, help = 'minimum number of cells a gene must be present in for it to be included as a potential highly variable gene')
	parser.add_argument('-preselectedGenesFile', default = None, type = str, help = 'path to a file containing a list of genes to use for PCA/clustering, expects genes to be entered as one per line with nothing else, set to None to skip')
	parser.add_argument('-dispersionMinMaxThreshold', nargs = 3, default = [0.025,4.0,0.5], type = float, metavar = ('min dispersion','max dispersion','dispersion threshold'), help = 'minimum, maximum and cutoff value for selecting highly dispersed genes')
	parser.add_argument('-bootstrapFrac', default = 0.8, type = float, help = 'fraction of cells to keep during bootstrapping')
	parser.add_argument('-regressOut', nargs = '+', default = ['n_counts','percent_mito'], type = str, metavar = 'first variable to regress out', help = 'variables to regress out from counts matrix')
	parser.add_argument('-usePCA', default = True, type = bool, help = 'whether to perform PCA prior to clustering')
	parser.add_argument('-kValue', default = 12, type = int, help = 'k value to use for constructing the nearest neighbor graph')
	parser.add_argument('-resolution', default = 1.0, type = float, help = 'resolution to use for calculating modularity')
	parser.add_argument('-clusterMin', default = 10, type = str, help = 'minimum number of cells that must be present in a cluster for it to be included in the final output')
	parser.add_argument('-verbose', default = True, type = bool, help = 'whether to print messages during processing')
	parser.add_argument('-trackIterations', default = True, type = bool, help = 'whether to track intermediate clustering results')
	parser.add_argument('-fileNameIteration', default = None, type = str, help = 'variable to allow appending numbers to file names for bootstrapping iterations')
	args = parser.parse_args()
	return args

def bootstrapClustering():
	args = parse_args()

	ex1 = experiment.Experiment(args.dataFile, args.outputLocation)
	ex1.filter(verbose = args.verbose, byBatch = args.byBatch, countsPercentileCutoffs = list(args.countsPercentileCutoffs),genesPercentileCutoffs = list(args.genesPercentileCutoffs), mitoPercentileMax = args.mitoPercentileMax, geneMin = args.geneMin)
	ex1.bootstrapCells(frac = args.bootstrapFrac)
	ex1.selectVariableGenes(preselectedGenesFile = args.preselectedGenesFile, dispersionMin = args.dispersionMinMaxThreshold[0], dispersionMax = args.dispersionMinMaxThreshold[1], dispersionThreshold = args.dispersionMinMaxThreshold[2])
	ex1.processData(regressOut=args.regressOut)
	if args.usePCA:
		ex1.selectPCs()
	ex1.computeNeighbors(kValue=args.kValue, usePCA=args.usePCA)
	ex1.cluster(resolution=args.resolution, clusterMin=args.clusterMin, trackIterations=args.trackIterations, fileNameIteration = args.fileNameIteration)

if __name__ == '__main__':
	bootstrapClustering()


