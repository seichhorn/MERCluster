import argparse
from MERCluster.classes import experiment
from MERCluster.classes import cluster_analysis

def parse_args():
	parser = argparse.ArgumentParser(description = 'process options for cluster analysis')

	parser.add_argument('clusterDir', type = str, help = 'path to directory containing results of full and bootstrap clustering')
	parser.add_argument('outputLocation', type = str, help = 'directory to write results')
	parser.add_argument('-plot', default = False, type = bool, help = 'whether to plot results in addition to writing result files' )
	parser.add_argument('-experimentData', type = str, help = 'path to raw data for constructing an experiment object')
	parser.add_argument('-cutUnstable', default = False, type = bool, help = 'remove cells assigned to unstable clusters')
	parser.add_argument('-geneIdentityFile', default = 'None', type = str, help = 'path to file containing genes to identify cell types')
	parser.add_argument('-cellType', default = None, type = str, help = 'name of cell type to perform analysis with, must match name in cluster files')
	args = parser.parse_args()

	return args


def processClustering():
	args = parse_args()

	cl1 = cluster_analysis.ClusterAnalysis(args.clusterDir, args.outputLocation, cellType = args.cellType)
	cl1.selectK(plot = args.plot)
	
	ex1 = experiment.Experiment(args.experimentData, args.outputLocation)
	cl1.identifyCellTypes(ex1, geneIdentityFile = args.geneIdentityFile, cutUnstable = args.cutUnstable, plot = args.plot)


if __name__ == '__main__':
	processClustering()


