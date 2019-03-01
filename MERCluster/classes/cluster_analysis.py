import os
from os.path import isfile, join

import sys
import numpy as np
import pandas as pd
import scanpy.api as sc
import datetime
import re
from collections import defaultdict
from sklearn import preprocessing 
import matplotlib as mpl
import matplotlib.pyplot as plt

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=150)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()


class ClusterAnalysis:

	def __init__(self, clusteringDir, outputDir, cellType = None):
		self.outputDir = outputDir
		if not os.path.exists(self.outputDir):
			os.makedirs(self.outputDir)

		if not os.path.exists(clusteringDir):
			print('The full clustering directory does not exist')
		else:
			self.clusteringDir = clusteringDir
		self.cellType = cellType


	def selectK(self, plot = True):

		kValuePattern = re.compile('kValue_\d+')
		resolutionPattern = re.compile('resolution_\d+')
		bootstrapPattern = re.compile('bootstrap_\d+')
		cellTypePattern = re.compile('type_[a-zA-Z]+')
		geneSetPattern = re.compile('geneset_[a-zA-Z0-9]+')

		baseDir = os.listdir(self.clusteringDir)
		onlyfiles = [f for f in baseDir if isfile(join(self.clusteringDir, f))]

		if self.cellType:
			cellTypes = ['type_'.format(self.cellType)]

		else:
			cellTypes = list(set([cellTypePattern.search(x).group() for x in onlyfiles if cellTypePattern.search(x)]))

		for cellType in cellTypes:
			filesToUse = [f for f in onlyfiles if cellType in f]

			geneSets = list(set([geneSetPattern.search(x).group() for x in onlyfiles if geneSetPattern.search(x)]))
			for geneSet in geneSets:
				filesToUseSubset = [f for f in filesToUse if geneSet in f]



				fullClustering = [x for x in filesToUseSubset if 'bootstrap' not in x]
				fullClustering = [[int(kValuePattern.search(x).group().split('_')[1]),int(resolutionPattern.search(x).group().split('_')[1]),x] for x in fullClustering]
				fullClustering.sort(key = lambda x: (x[0],x[1]))
				fullClustering = [x[2] for x in fullClustering]

				for f in fullClustering:
					if f.split('.')[-1] == 'csv':
						currentFile = pd.read_csv(self.clusteringDir+f,index_col = 0)
					elif f.split('.')[-1] == 'txt':
						currentFile = pd.read_table(self.clusteringDir+f,index_col = 0)
					else:
						print('unknown file extension {}'.format(f.split('.')[-1]))
					kval = kValuePattern.search(f).group().split('_')[1]
					res = resolutionPattern.search(f).group().split('_')[1]
					currentFile.columns = ['{}_{}'.format(kval,res)]
					if f == fullClustering[0]:
						fullClusteringDF = currentFile.copy(deep=True)
					else:
						fullClusteringDF = fullClusteringDF.merge(currentFile,left_index = True, right_index = True)
						

				bootstrapClustering = [x for x in filesToUseSubset if 'bootstrap' in x]
				bootstrapClustering = [[int(bootstrapPattern.search(x).group().split('_')[1]),int(kValuePattern.search(x).group().split('_')[1]),int(resolutionPattern.search(x).group().split('_')[1]),x] for x in bootstrapClustering]
				bootstrapClustering.sort(key = lambda x: (x[0],x[1],x[2]))
				bootstrapClustering = [x[3] for x in bootstrapClustering]


				kValues = fullClusteringDF.columns.values.tolist()
				results = []
				stableClusters = dict()
				for kValue in kValues:   
					currentFull = fullClusteringDF.loc[:,[kValue]]
					totalLen = currentFull.shape[0]
					currentFull = currentFull[currentFull[kValue] != -1]
					currentFull.columns = ['Full clustering']
					bootstrapCollection = []
					hit = 0
					specificK = re.compile('kValue_{}_resolution_{}_'.format(kValue.split('_')[0],kValue.split('_')[1]))
					for f in bootstrapClustering:
						if specificK.search(f):
							if f.split('.')[-1] == 'csv':
								currentFile = pd.read_csv(self.clusteringDir+f,index_col = 0)
							elif f.split('.')[-1] == 'txt':
								currentFile = pd.read_table(self.clusteringDir+f,index_col = 0)
							else:
								print('unknown file extension {}'.format(f.split('.')[-1]))
							bootstrapCollection.append(currentFile)
					for boot in bootstrapCollection:
						currentBoot = boot
						currentBoot = currentBoot[currentBoot['kValue_{}_resolution_{}'.format(kValue.split('_')[0],kValue.split('_')[1])] != -1]
						currentBoot.columns = ['Bootstrap clustering']
						combinedFullBoot = currentFull.merge(currentBoot,left_index = True, right_index = True)

						recovery = combinedFullBoot.groupby(['Full clustering','Bootstrap clustering']).size().unstack().max(1).div(combinedFullBoot.groupby('Full clustering').size())
						if hit == 0:
							recoveryDF = pd.DataFrame(recovery)
							hit+=1
						else:
							recoveryDF = pd.concat([recoveryDF,pd.DataFrame(recovery)],axis = 1)
					stableEntries = len(recoveryDF[recoveryDF.median(1) > 0.5])
					totalEntries = len(recoveryDF)
					recoveredCells = len(currentFull[currentFull['Full clustering'].isin(recoveryDF[recoveryDF.median(1) >= 0.5].index.values.tolist())])/totalLen
					stableClusters[kValue] = recoveryDF[recoveryDF.median(1) >= 0.5].index.values.tolist()
					results.append([kValue,stableEntries,totalEntries,recoveredCells])

					clusteringResult = fullClusteringDF.loc[:,[kValue]].copy()
					clusteringResult['Stable'] = clusteringResult[kValue].isin(stableClusters[kValue]) 

					if not os.path.exists(str(self.outputDir) + 'all_analyses/'):
						os.makedirs(str(self.outputDir) + 'all_analyses/')

					clusteringResult.to_csv(str(self.outputDir)+'all_analyses/stability_analysis_kValue_{}_resolution_{}_type_{}_geneset_{}.txt'.format(kValue.split('_')[0],kValue.split('_')[1],cellType.split('_')[1],geneSet.split('_')[1]),sep = '\t')


				position = 0
				secondCheck = []
				while position < len(results):
					if results[position][-1] >= 0.9:
						secondCheck.append(position)
					position += 1

				toCheck = [results[x] for x in secondCheck]
				toCheck.sort(key = lambda x: x[1], reverse = True)
				selectedK = str(toCheck[0][0])

				print('Choosing {}'.format(selectedK))
				print(results)

				clusteringResult = fullClusteringDF.loc[:,[selectedK]].copy()
				clusteringResult['Stable'] = clusteringResult[selectedK].isin(stableClusters[selectedK]) 

				clusteringResult.to_csv(str(self.outputDir)+'selected_stability_analysis_kValue_{}_resolution_{}_type_{}_geneset_{}.txt'.format(selectedK.split('_')[0],selectedK.split('_')[1],cellType.split('_')[1],geneSet.split('_')[1]),sep = '\t')


				self.clusterLabels = clusteringResult.iloc[:,[0]]
				self.stability = clusteringResult.iloc[:,[1]]

				if plot:
					f,ax = plt.subplots(1,1,figsize = (10,10))
					plt.scatter([x[3] for x in results],[x[1] for x in results], cmap = 'tab20')
					for i, txt in enumerate([x[0] for x in results]):
					    ax.annotate(str(txt), ([x[3] for x in results][i],[x[1] for x in results][i]))
					ax.set_xlim(left = 0, right = 1)
					f.savefig(str(self.outputDir)+'selected_stability_analysis_type_{}_geneset_{}.png'.format(cellType,geneSet),bbox_inches = 'tight')


	def identifyCellTypes(self, experiment, geneIdentityFile = 'None', cutUnstable = False, plot = True):	
		# The order of cell types in the geneIdentityFile determines the hierarchy for breaking ties if a cell is enriched
		# for more than one marker type

		# If a file for identifying gene identity isn't supplied it uses the following by default

		# Excitatory neuronal: Slc17a6, Slc17a7
		# Inhibitory neuronal: Gad1, Gad2, Slc32a1
		# Astrocytic: Aqp4, Mlc1
		# Microglial: Selplg, Slc15a3
		# Endothelial: Fn1, Slco1a4
		# Pericyte: Myh11, Lmod1
		# Pan oligodendrocyte: Gjc3, Sgk1

		if geneIdentityFile == 'None':
			genesForIdentification = defaultdict(list)
			genesForIdentification['Inhibitory_neuron'] = ['Gad1','Slc32a1']
			genesForIdentification['Excitatory_neuron'] = ['Slc17a6','Slc17a7']
			genesForIdentification['Oligodendrocyte'] = ['Gjc3','Sgk1']
			genesForIdentification['Astrocyte'] = ['Aqp4','Mlc1']
			genesForIdentification['Microglia'] = ['Selplg','Slc15a3']
			genesForIdentification['Endothelial_cell'] = ['Fn1','Slco1a4']
			genesForIdentification['Pericyte'] = ['Myh11','Lmod1'] 

			cellTypes = ['Inhibitory_neuron','Excitatory_neuron','Oligodendrocyte','Astrocyte','Microglia','Endothelial_cell','Pericyte']
			allGenes = [x for y in list(genesForIdentification.values()) for x in y]
		else:
			cellTypes = []
			allGenes = []
			genesForIdentification = defaultdict(list)
			with open(geneIdentityFile,'r') as geneIdentityOpen:
				for line in geneIdentityOpen:
					cellType,gene = line.rstrip('\n').split('\t')
					genesForIdentification[cellType].append(gene)
					if cellType not in cellTypes:
						cellTypes.append(cellType)
					allGenes.append(gene)

#TODO: add a check to make sure that the genes you're requesting are in the data to begin with, would be helpful to print the missing ones
		rawData = pd.DataFrame(data = experiment.dataset[:,allGenes].X, index = experiment.dataset.obs.index, columns = allGenes)
		
		clFiles = os.listdir(self.clusteringDir)
		finalFiles = [f for f in clFiles if 'selected_stability_analysis' in f]

		cellTypePattern = re.compile('type_[a-zA-Z]+')
		geneSetPattern = re.compile('geneset_[a-zA-Z0-9]+')
		
		geneSets = list(set([geneSetPattern.search(x).group() for x in finalFiles if geneSetPattern.search(x)]))
		for geneSet in geneSets:
			filesToUseSubset = [f for f in finalFiles if geneSet in f]

			for f in filesToUseSubset:
				cellType = cellTypePattern.search(f).group().split('_')[-1]
				temp = pd.read_table(str(self.clusteringDir)+str(resultFile),index_col = 0)
				temp.iloc[:,0] = str(cellType) + temp.iloc[:,0].astype(str)
				if f == filesToUseSubset[0]:
					clResults = temp.copy(deep=True)
				else:
					clResults = pd.concat([clResults,temp])



			clusterData = rawData.merge(clResults.iloc[:,[0]], left_index = True, right_index = True)
			if clusterData.shape[0] < rawData.shape[0]:
				print('Some cells were dropped, likely because they were absent from the cluster labels')
			if clusterData.shape[0] > rawData.shape[0]:
				print('There are likely redundant index identifiers in the dataset, this should be fixed')

			if cutUnstable:
				clusterData = clusterData.merge(clResults.iloc[:,[1]], left_index = True, right_index = True)
				clusterData = clusterData[clusterData['Stable']].iloc[:,:-1]

			clusterMeans = clusterData.groupby(clusterData.columns.values.tolist()[-1]).mean()
			clusterZScores = pd.DataFrame(data = preprocessing.scale(clusterMeans), index = clusterMeans.index, columns = clusterMeans.columns)

			accounted = []

			with open(str(self.outputDir) + 'cellTypes_analysis_geneset_{}.txt'.format(geneSet),'w') as outPathOpen:
				for cellType in cellTypes:
					genes = genesForIdentification[cellType]
					mean = clusterZScores.loc[:,genes]
					hits = [x for x in mean[mean>0.0].dropna(how = 'all',axis = 0).index.values.tolist() if x not in accounted]
					accounted.extend(hits)
					for x in hits:
						outPathOpen.write('{}\t{}\n'.format(x, cellType))

				unknown = [x for x in clusterZScores.index.values.tolist() if x not in accounted]
				for x in unknown:
					outPathOpen.write('{}\t{}\n'.format(x, 'unknown'))			

			if plot:
				f,axs = plt.subplots(1,1,figsize = (20,20))
				sns.heatmap(clusterZScores,vmax = 1, vmin = 0, cmap = 'Purples')
				f.savefig(tr(self.outputDir) + 'cellTypes_analysis_geneset_{}.png'.format(geneSet),bbox_inches = 'tight', dpi = 200)


