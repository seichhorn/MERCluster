# MERCluster

Code to simplify clustering of sc-seq and MERFISH data, largely based on scanpy functionality


Install:  
Code was written based on the development version 1.3 of scanpy, this originally offered additional functionality to that available from conda. Functions have not yet been updated to accomodate the current scanpy release. To recapitulate the code base I'm using:

conda create -n scanpy  
source activate scanpy  

clone MERCluster  
cd MERCluster/  
pip install -e .  
  
conda install seaborn scikit-learn statsmodels numba  
git clone git@github.com:theislab/scanpy.git  
git checkout 37851434b2077be2debcb2c4ffd294f6d0ac9707  
cd scanpy  
pip install -e .  
conda install -c conda-forge python-igraph louvain  
conda install notebook ipykernel  
ipython kernel install --user --name scanpy   
pip install dca  
pip install -U tensorflow  
git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git  
cd Multicore-TSNE/  
pip install .  
  

When I load scanpy and ask the versions I get:  
scanpy==1.3.2+85.g3785143 anndata==0.6.12 numpy==1.14.5 scipy==1.1.0 pandas==0.23.4 scikit-learn==0.19.1 statsmodels==0.9.0 python-igraph==0.7.1 louvain==0.6.1 
