# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 21:28:45 2017

@author: ashutosh
"""

import os
import sys
import pickle
import command_args
import create_network
import numpy as np
import matplotlib.pyplot as plt

#xtc,tpr,ndx,skip,out=command_args.parseargs()
xtc,tpr,ndx,skip,b,e,out=command_args.parseargs()
if os.path.isfile("all_contacts.txt"):
    print("Contacts already existing in all_contacts.txt")
    with open('all_contacts.txt','rb') as all_con:
        all_contacts=pickle.load(all_con)
if os.path.isfile("contact_matrix.txt"):
    print("Contact matrix already exixts in folder. Now reading it")
    contact_mat=np.loadtxt('contact_matrix.txt')
print(all_contacts)
    
#else:
#    contact_mat=np.zeros([len(all_contacts),5])
#    count=0
#    for i1 in range(0,1000,int(skip)):
#        with open('frame'+str(i1)+'contacts.txt','rb') as rcon:
#            con_read=pickle.load(rcon)
#        for item in con_read:
#            if item in list(all_contacts):
#                index=list(all_contacts).index(item)
#                contact_mat[index,count]=1
#        count=count+1
#    np.savetxt('contact_matrix.txt',contact_mat,fmt="%d")
#    
#plt.imshow(contact_mat,aspect="auto",interpolation='nearest')
#[num_cont,num_frames]=np.shape(contact_mat)
#sum_contact=np.sum(contact_mat,axis=1)
#cont_prob=sum_contact/num_frames    #Calculate contact probability
#plt.plot(range(len(cont_prob)),sorted(cont_prob,reverse=True))
#plt.show()         