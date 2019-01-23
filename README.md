# MERCluster

Code to simplify clustering of sc-seq and MERFISH data, largely based on scanpy functionality


Install:
Code was written based on the development version of scanpy, this originally offered additional functionality to that available from conda. To recapitulate the code base I'm using:

clone MERClusters  
cd MERCluster/  
conda create --name myenv --file /path/to/spec-file.txt  
source activate myenv  

clone scanpy from https://github.com/theislab/scanpy  
cd scanpy/  
pip install -e .  

git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git  
cd Multicore-TSNE/  
pip install .  

