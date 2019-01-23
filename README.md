# MERCluster

Code to simplify clustering of sc-seq and MERFISH data, largely based on scanpy functionality


Install:  
Code was written based on the development version of scanpy, this originally offered additional functionality to that available from conda. To recapitulate the code base I'm using:

clone MERCluster  
cd MERCluster/  
conda create --name myenv --file /path/to/spec-file.txt  
source activate myenv  

clone scanpy from https://github.com/theislab/scanpy  
cd scanpy/  
pip install -e .  

git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git  
cd Multicore-TSNE/  
pip install .  

When I load scanpy and ask the versions I get:  
scanpy==1.3.2+85.g3785143 anndata==0.6.12 numpy==1.14.5 scipy==1.1.0 pandas==0.23.4 scikit-learn==0.19.1 statsmodels==0.9.0 python-igraph==0.7.1 louvain==0.6.1 
