# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 10:56:46 2017
This script calculates the dynamic network from trajectory 
@author: ashutosh
"""

import os
import subprocess
import shutil
from distutils import spawn
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser # for parsing PDB file
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale

#import command_args
def extract_frame(xtc,tpr,ndx,i):
    print("Trajectory file read: %s" %(xtc))
    print("Structure file read: %s" %(tpr))
    print("Index file read: %s" %(ndx))    
    choice=open('choice.txt','r')                            # File with option to write the frames
    fr=str(i)                                                #Frame number to be dumped converted to string for use within subprocess
    fname='frame'+fr+'.pdb'                                  #PDB file in which frame is dumped
    sout=open(os.devnull,'w')                                #dumping standard output from gromacs
    if spawn.find_executable("gmx"):
        subprocess.call(['gmx','trjconv', '-f', xtc, '-s', tpr, '-dump', fr,'-o', fname, '-n', ndx],stdin=choice, stdout=sout,stderr=subprocess.STDOUT)
    else:
        subprocess.call(['trjconv', '-f', xtc, '-s', tpr, '-dump', fr,'-o', fname, '-n', ndx],stdin=choice,stderr=subprocess.STDOUT)
    choice.close()
    print("Frame %s dumped in file %s" %(fr,fname))
    return fname

def find_contact(fname):
    if os.path.isfile(fname)==True:                                    #checks if the name saved in file_list is a file
        print(fname)
        parser=PDBParser(PERMISSIVE=1)                                     #parser for PDB file
        temp_struct=parser.get_structure(fname[0:-4], fname)                # parsing of PDB file
        model=temp_struct[0]
        contact_count=0
        contacts=[]
        resnum=0
        res_ids=[]
        chain_seq=[]
        for res in model.get_residues():
            if PDB.is_aa(res):
                res_ids.append(res.get_id()[1])
                chain_seq.append(res.get_parent().get_id())
                resnum=resnum+1
        for i in range(0,len(res_ids)):
            for j in range(0,len(res_ids)):
                atom1=temp_struct[0][chain_seq[i]][res_ids[i]]['CA']
                atom2=temp_struct[0][chain_seq[j]][res_ids[j]]['CA']
                if (atom1-atom2 <= 7) and (abs(int(res_ids[i])-res_ids[j])>2):
                    contact=(res_ids[i],chain_seq[i],res_ids[j],chain_seq[j])
                    contact_rev=(res_ids[j],chain_seq[j],res_ids[i],chain_seq[i])
                    if contact in contacts or contact_rev in contacts :
                        pass
                    else:
                        #contacts.append((atom1.get_id(),resid1.get_id()[1],atom2.get_id(),resid2.get_id()[1]))
                        contacts.append(contact)
                        contact_count=contact_count+1
                else:
                    pass 
    return contacts 
                
            
            

# def find_contact(fname):
#     if os.path.isfile(fname)==True:                                    #checks if the name saved in file_list is a file
#         print(fname)
#         parser=PDBParser(PERMISSIVE=1)                                     #parser for PDB file
#         temp_struct=parser.get_structure(fname[0:-4], fname)                # parsing of PDB file
#         model=temp_struct[0]
#         contact_count=0
#         contacts=[]
#         for residue1 in model.get_residues():                      #loop over all the residues in model    
#             chain1=residue1.get_parent()                           #chain id for first loop
#             if chain1.get_id()=='A' or chain1.get_id()=='B':       # residues belonging only to chain A or B processed
#                 for residue2 in model.get_residues():
#                     chain2=residue2.get_parent()  
#                     if chain2.get_id()=='A' or chain2.get_id()=='B':
#                         temp_resid1=residue1.get_id()
#                         temp_resid2=residue2.get_id()
#                         if temp_resid1[0]==" " and temp_resid2[0]==" ": # only aa residues, excludes Water and other het atoms
#                             if residue1['CA'] and residue2['CA']:       #selection of CA atoms
#                                 atom1=residue1['CA']
#                                 atom2=residue2['CA']
#                                 if (atom1-atom2 <= 7) and (abs(int(residue1.get_id()[1])-residue2.get_id()[1])>2):
#                                     contact=(residue1.get_id()[1],chain1.get_id(),residue2.get_id()[1],chain2.get_id())
#                                     contact_rev=(residue2.get_id()[1],chain2.get_id(),residue1.get_id()[1],chain1.get_id())
#                                     if contact in contacts or contact_rev in contacts :
#                                         pass
#                                     else:
#                                             #contacts.append((atom1.get_id(),resid1.get_id()[1],atom2.get_id(),resid2.get_id()[1]))
#                                         contacts.append(contact)
#                                         contact_count=contact_count+1
#                                 else:
#                                     pass 
#                         else:
#                             pass
#     return contacts 
    

def extract_dynamic_contacts(low_cut,high_cut,all_contacts,contact_mat,num_frames):
    sum_contact=np.sum(contact_mat,axis=1)
    cont_prob=sum_contact/num_frames    #Calculate contact probability
    contacts_selected=[]
    selected_cont_rows=[]
    ind=0
    for i in cont_prob:
        if low_cut<i<high_cut:
            contacts_selected.append(all_contacts[ind])
            selected_cont_rows.append(contact_mat[ind,:])
            #print(i,ind)
        ind=ind+1
    selected_cont_mat=np.vstack(selected_cont_rows)
    return contacts_selected,selected_cont_mat

def PCA_contact_mat(selected_cont_mat):

#Load data set
#data = pd.read_csv('Big_Mart_PCA.csv')

#convert it to numpy arrays
#X=data.values

#Scaling the values
#X = scale(X)

    pca = PCA(n_components=2)
    pca.fit(selected_cont_mat)
    trans_cont_mat=pca.transform(selected_cont_mat)
    
    #The amount of variance that each PC explains
    #var= pca.explained_variance_ratio_
    
    #Cumulative Variance explains
    var1=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
    print(var1)
   # plt.plot(var1)
    plt.scatter(trans_cont_mat[:,0],trans_cont_mat[:,1])
    plt.show()
        
        
        
    