# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 12:59:53 2017

@author: ashutosh
"""
import command_args
#import numpy as np
#import matplotlib.pyplot as plt
#
#contact_mat=np.loadtxt('contact_matrix.txt')
#sum_contact=np.sum(contact_mat,axis=1)
##print(np.shape(sum_contact))
##print
#plt.bar(range(len(sum_contact)),sorted(sum_contact,reverse=True))
#
##plt.imshow(contact_mat,aspect="auto",interpolation='nearest')
#plt.show()
xtc,tpr,ndx,skip,b,e,out=command_args.parseargs()
for i1 in range(int(b),int(e),int(skip)):
        #with open('frame'+str(i1)+'contacts.txt','rb') as rcon:
            #con_read=pickle.load(rcon)
#        with open('frame'+str(i1)+'contacts.json','r') as rcon:
#            con_read=json.loads(rcon)
        with open('frame'+str(i1)+'contacts.txt','r') as rcon:
            con_read=[tuple(i_c.strip().split(' ')) for i_c in rcon]
            print(con_read)