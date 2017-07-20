'''
Created on 08-Jun-2015
Calculates CA based distance matrices for the PDB files in folder
@author: ashutosh
'''
import os
import shutil
import numpy as np
from Bio.PDB.PDBParser import PDBParser # for parsing PDB file
import adjacency_mat as admat
import residue_label_extract as rle
from numpy import sqrt
def distance_calc(p,fold):
    file_list=os.listdir(fold)
    #print file_list
    os.chdir(fold)
    parser=PDBParser(PERMISSIVE=1)                                     #parser for PDB file
    for f in file_list:
        if os.path.isfile(f)==True:                                    #checks if the name saved in file_list is a file
            print f
            temp_struct=parser.get_structure(f[0:4], f)                # parsing of PDB file
            model=temp_struct[0]                                       # parsing only first model
            residue_labels=rle.residue_label_extract(model)
            #count1=0
            dist_at1_at2=[]                                             #intializing distance array
            #count3=0
            for residue1 in model.get_residues():                      #loop over all the residues in model    
                chain1=residue1.get_parent()                           #chain id for first loop
                if chain1.get_id()=='A' or chain1.get_id()=='B':       # residues belonging only to chain A or B processed
                    #count1=count1+1
                    #print count1
                    #count2=0
                    for residue2 in model.get_residues():
                        chain2=residue2.get_parent()  
                        if chain2.get_id()=='A' or chain2.get_id()=='B':
                            #count2=count2+1
                            
                            #print count2
                            temp_resid1=residue1.get_id()
                            temp_resid2=residue2.get_id()
                            if temp_resid1[0]==" " and temp_resid2[0]==" ": # only aa residues, excludes Water and other het atoms
                                if residue1['CA'] and residue2['CA']:       #selection of CA atoms
                                    #count3=count3+1
                                    atom1=residue1['CA']
                                    atom2=residue2['CA']
                                    dist_at1_at2.append(atom1-atom2)        # distance calculation and addition to array
                                else:
                                    pass 
                            else:
                                pass
                                              
        #print dist_at1_at2
        #print count3
        print sqrt(len(dist_at1_at2)) 
         
        dist_mat=np.reshape(dist_at1_at2,(sqrt(len(dist_at1_at2)),sqrt(len(dist_at1_at2))))
        os.mkdir(f[0:4])
        temp_f=f[0:4]#+'/'+f
        shutil.move(f,temp_f)
        temp_dist_file='./'+f[0:4]+'/'+f[0:4]+'_distmat.txt'
        temp_label_file='./'+f[0:4]+'/'+'labels_'+f[0:4]+'.txt'
        np.savetxt(temp_dist_file, dist_mat, fmt="%f",delimiter='\t')
        np.savetxt(temp_label_file,residue_labels,fmt="%d",delimiter='\t')
        temp_adj_file='./'+f[0:4]+'/'+f[0:4]+'_adjacency.txt'
        admat.adjacency_mat(dist_mat,temp_adj_file,7)  
        
        temp_struct=None
        dist_mat=None
        temp_dist_file=None
        temp_adj_file=None
        #----------------------------------------------------------- count1=None
        #----------------------------------------------------------- count2=None
        #----------------------------------------------------------- count3=None
        chain1=None
        chain2=None
        residue1=None
        residue2=None
        temp_resid1=None
        temp_resid2=None   
        model=None
        dist_at1_at2=None
        atom1=None
        atom2=None   
    return