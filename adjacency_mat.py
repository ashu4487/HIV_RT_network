'''
Created on 09-Jun-2015
Creates adjacency matrix from distance matrix using cut-off
@author: ashutosh
'''
import numpy as np
def adjacency_mat(dist_mat,file_name,threshold):
    mat_size=dist_mat.shape
    adj=np.zeros([mat_size[0],mat_size[1]])
    for i in range(0,mat_size[0]):
        for j in range(0,mat_size[1]):
            if dist_mat[i][j]<threshold and i != j:
                adj[i][j]=1
            else:
                adj[i][j]=0
    np.savetxt(file_name,adj,fmt="%d",delimiter='\t')
    mat_size=None
    return