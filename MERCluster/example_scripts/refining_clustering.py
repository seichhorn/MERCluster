import numpy as np
import pandas as pd
import seaborn as sns
import os
import matplotlib as mpl
import scanpy.api as sc
import datetime
from sklearn import preprocessing 
import sys
sys.path.append('/n/home13/seichhorn/Python_code/MERCluster/')
from classes import experiment
from classes import cluster_analysis
from scipy.spatial import cKDTree
import scipy 
import networkx 
from networkx.algorithms.components.connected import connected_components

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=150)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()
sc.settings.figdir = '/n/home13/seichhorn/Hypothalamus/'

today = datetime.date.today()

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

cellType = sys.argv[1]

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    


rawData = '/n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/raw_data/combined_data_geneNames_180912.h5ad'
outPath = '/n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/190125/Round2/clustering_revisions/'
clusterPath = '/n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/190125/Round2/clustering/'

ex1 = experiment.Experiment(rawData,outPath)

ex1.cutToCellList('/n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/190125/Round1/analysis/cellTypes.txt','/n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/190125/Round1/analysis/stability_analysis.txt', cellType)

ex1.selectVariableGenes(preselectedGenesFile = 'None', dispersionMin = 0.025, dispersionMax = 4.0, dispersionThreshold = 0.5)

regressOut=['n_counts','percent_mito']
sc.pp.log1p(ex1.dataset)

if 'n_counts' in regressOut:
        if 'n_counts' not in ex1.dataset.obs.columns.values.tolist():
            self.dataset.obs['n_counts'] = ex1.dataset.X.sum(1)
if 'percent_mito' in regressOut:
    if 'percent_mito' not in ex1.dataset.obs.columns.values.tolist():
        mito_genes = [name for name in ex1.dataset.var_names if name.startswith('mt-')]
        ex1.dataset.obs['percent_mito'] = np.sum(ex1.dataset[:, mito_genes].X, axis=1) / np.sum(ex1.dataset.X, axis=1)

sc.pp.regress_out(ex1.dataset, regressOut)


cl1 = cluster_analysis.ClusterAnalysis(clusterPath, outPath, cellType = cellType)

from os.path import isfile, join

baseDir = os.listdir(cl1.clusteringDir)
onlyfiles = [f for f in baseDir if isfile(join(cl1.clusteringDir, f))]

if cl1.cellType:
    onlyfiles = [f for f in onlyfiles if 'type_{}'.format(cl1.cellType) in f]


for f in onlyfiles:
    pathToClResults = cl1.clusteringDir+f
    clResults = pd.read_table(pathToClResults,index_col = 0)

    merged = clResults.merge(pd.DataFrame(ex1.dataset.X, index = ex1.dataset.obs.index, columns = ex1.dataset.var.index), left_index = True, right_index= True)
    means = merged.groupby(merged.columns.values.tolist()[0]).mean()

    indexes = means.index.values.tolist()
    indexes = [x for x in indexes if x != -1]
    starting = len(indexes)
    tree = cKDTree(means.loc[indexes,:])
    a,b = tree.query(means.loc[indexes,:],k=2)
    pairs = [list(x) for x in b]

    tracking = []
    toMerge = []
    for pair in pairs:
        first,second = pair
        g1 = merged[merged[merged.columns.values.tolist()[0]] == first].iloc[:,1:]
        g2 = merged[merged[merged.columns.values.tolist()[0]] == second].iloc[:,1:]
        [t,p] = scipy.stats.ttest_ind(g1,g2)
        fc = abs(np.log2(abs(g1.mean()/g2.mean())))
        pDF = pd.DataFrame(p).reset_index()
        pDF = pDF.sort_values(by = 0).reset_index(drop = True)

        i = list(range(1,pDF.shape[0]+1))
        pDF['adjP'] = [(pDF.shape[0]-i[x]+1)*pDF.loc[x,0] for x in range(len(p))]
        pDF['passing_pVal'] = pDF['adjP']<=0.01/len(pairs)
        pDF['fc'] = fc[pDF['index'].values.tolist()].values.tolist()
        pDF['passing_fc'] = pDF['fc'] >= 1

        kept = pDF[pDF['passing_pVal'] & pDF['passing_fc']].copy(deep=True)
        if len(kept) == 0:
            toMerge.append(pair)
        else:
            mask = kept.adjP < 1*(10**-20)
            column_name = 'adjP'
            kept.loc[mask, column_name] = 1*(10**-20)
            kept['invertedP'] = -np.log10(kept['adjP'])
            toBeat = kept.sum()['invertedP']
            tracking.append(toBeat)
            if toBeat<150:
                toMerge.append(pair)

    if len(toMerge) > 0:
        while len(toMerge) > 0:
            G = to_graph(toMerge)
            toMergeConnected = sorted(connected_components(G), key = len, reverse=True)
            toMergeConnected = [list(x) for x in toMergeConnected]
            presentInMerge = [x for y in toMergeConnected for x in y]
            allSets = toMergeConnected + [[x] for x in indexes if x not in presentInMerge]
            namingDict = {}
            counter = 0
            for element in allSets:
                for e in element:
                    namingDict[e] = counter
                counter += 1
            namingDict[-1] = -1
            mergedUpdate = merged.copy(deep=True)
            mergedUpdate['recast'] = mergedUpdate[mergedUpdate.columns.values.tolist()[0]].map(namingDict)
            mergedUpdate['recast'] = mergedUpdate['recast'].values.astype(int)

            mergedUpdate = mergedUpdate.iloc[:,[-1]+list(range(1,mergedUpdate.shape[1]-1))]
            mergedUpdate.columns = ['updated_kValues'] + mergedUpdate.columns.values.tolist()[1:]    

            means = mergedUpdate.groupby(mergedUpdate.columns.values.tolist()[0]).mean()

            indexes = means.index.values.tolist()
            indexes = [x for x in indexes if x != -1]

            tree = cKDTree(means.loc[indexes,:])

            a,b = tree.query(means.loc[indexes,:],k=2)
            pairs = [list(x) for x in b]

            toMerge = []
            for pair in pairs:
                first,second = pair
                g1 = mergedUpdate[mergedUpdate[mergedUpdate.columns.values.tolist()[0]] == first].iloc[:,1:]
                g2 = mergedUpdate[mergedUpdate[mergedUpdate.columns.values.tolist()[0]] == second].iloc[:,1:]
                [t,p] = scipy.stats.ttest_ind(g1,g2)
                fc = abs(np.log2(abs(g1.mean()/g2.mean())))
                pDF = pd.DataFrame(p).reset_index()
                pDF = pDF.sort_values(by = 0).reset_index(drop = True)

                i = list(range(1,pDF.shape[0]+1))
                pDF['adjP'] = [(pDF.shape[0]-i[x]+1)*pDF.loc[x,0] for x in range(len(p))]
                pDF['passing_pVal'] = pDF['adjP']<=0.01/len(pairs)
                pDF['fc'] = fc[pDF['index'].values.tolist()].values.tolist()
                pDF['passing_fc'] = pDF['fc'] >= 1

                kept = pDF[pDF['passing_pVal'] & pDF['passing_fc']].copy(deep=True)

                if len(kept) == 0:
                    toMerge.append(pair)
                else:
                    mask = kept.adjP < 1*(10**-20)
                    column_name = 'adjP'
                    kept.loc[mask, column_name] = 1*(10**-20)
                    kept['invertedP'] = -np.log10(kept['adjP'])
                    toBeat = kept.sum()['invertedP']
                    tracking.append(toBeat)
                    if toBeat<150:
                        toMerge.append(pair)

            merged = mergedUpdate.copy(deep=True)

        sorted_index =  merged.groupby('updated_kValues').size().sort_values(ascending = False).index.values.tolist()
        sorted_index = [x for x in sorted_index if x != -1]
        ending = len(sorted_index)
        conversionDict = dict(zip(sorted_index,list(range(len(sorted_index)))))
        conversionDict[-1] = -1
        merged['sorted_kValues'] = merged['updated_kValues'].map(conversionDict)
        forOut = merged.loc[:,['sorted_kValues']]
        forOut.columns = clResults.columns.values.tolist()
        forOut.to_csv(outPath+f)
        print('started with {} ended with {}'.format(starting,ending))
    else:
        clResults.to_csv(outPath+f)
        print('started with {} ended with {}'.format(starting,starting))

