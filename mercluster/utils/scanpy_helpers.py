import numpy as np
import scipy.sparse as sp
import pandas as pd


def jaccard_kernel(sparseConnectivites):
	"""
	This function is directly copied from https://github.com/jacoblevine/PhenoGraph

	Compute Jaccard coefficient between nearest-neighbor sets
	"""
	n = sparseConnectivites.shape[0]
	s = list()
	r = list()
	j = list()
	for i in range(n):
		shared_neighbors = np.fromiter((len(
			set(sparseConnectivites[i].indices).intersection(
				set(sparseConnectivites[j].indices))) for j in
										sparseConnectivites[i].indices),
									   dtype=float)
		num_neighbors = np.fromiter((len(
			set(sparseConnectivites[i].indices)) + len(
			(set(sparseConnectivites[j].indices))) for j in
									 sparseConnectivites[i].indices),
									dtype=float)
		s.extend(shared_neighbors / (num_neighbors - shared_neighbors))
		r.extend([i] * len(sparseConnectivites[i].indices))
		j.extend(sparseConnectivites[i].indices)
	return r, j, s


def neighbor_graph(kernel, connectivities, directed=False, prune=False):
	"""
	This function is directly copied from https://github.com/jacoblevine/PhenoGraph

	Compute neighbor graph based on supplied kernel and connectivities
	"""

	r, j, s = kernel(connectivities)
	graph = sp.coo_matrix((s, (r, j)), shape=(
	connectivities.shape[0], connectivities.shape[0]))

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

def shuffler(array: np.array):
	"""
	returns a version of a numpy array where the values in each column
	have been shuffled independently
	Args:
		array: numpy array to use in shuffling
	Returns:
		a shuffled numpy array
	"""
	idx = [np.random.choice(array.shape[0], array.shape[0], replace=False) for
		   _ in range(array.shape[1])]

	holding = array[np.array(idx).T, np.arange(array.shape[1])]
	return holding
