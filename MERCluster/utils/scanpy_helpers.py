#Code directly ported over from phenograph code

import numpy as np
import scipy.sparse as sp
import pandas as pd

def jaccard_kernel(sparseConnectivites):
	"""
	Compute Jaccard coefficient between nearest-neighbor sets
	:param idx: numpy array of nearest-neighbor indices
	:return (i, j, s): tuple of indices and jaccard coefficients, suitable for constructing COO matrix
	"""
	n = sparseConnectivites.shape[0]
	s = list()
	r = list()
	j = list()
	for i in range(n):
		shared_neighbors = np.fromiter((len(set(sparseConnectivites[i].indices).intersection(set(sparseConnectivites[j].indices))) for j in sparseConnectivites[i].indices),dtype = float)
		num_neighbors = np.fromiter((len(set(sparseConnectivites[i].indices))+len((set(sparseConnectivites[j].indices))) for j in sparseConnectivites[i].indices),dtype = float)
		s.extend(shared_neighbors / (num_neighbors - shared_neighbors))
		r.extend([i]*len(sparseConnectivites[i].indices))
		j.extend(sparseConnectivites[i].indices)
	return r, j, s

def neighbor_graph(kernel,connectivities, directed = False, prune = False):
	r,j,s = kernel(connectivities)
	graph = sp.coo_matrix((s, (r, j)),shape = (connectivities.shape[0],connectivities.shape[0]))

	if not directed:
		if not prune:
			# symmetrize graph by averaging with transpose
			sg = (graph + graph.transpose()).multiply(.5)
		else:
			# symmetrize graph by multiplying with transpose
			sg = graph.multiply(graph2.transpose())
		# retain lower triangle (for efficiency)
		graph = sp.tril(sg, -1)
	
	return graph.tocsr()

def minimum_cluster_size(communities, min_size = 10):
	'''
	Takes a pandas dataframe with cells as the index, and cluster identifier as the column, and removes
	all clusters smaller than minimum size
	'''
	communities.columns = ['louvain']
	communitiesGrouped = communities.groupby('louvain').size()

	toZero = communitiesGrouped[communitiesGrouped<int(min_size)].index.values.tolist()
	mask = communities['louvain'].isin(toZero)
	communities['louvain'] = communities['louvain'].where(~mask,other = -1)

	return communities


def shuffler(matrix):
    idx = [np.random.choice(matrix.shape[0], matrix.shape[0], replace=False) for x in range(matrix.shape[1])]
    
    holding = matrix[np.array(idx).T, np.arange(matrix.shape[1])]
    return holding


