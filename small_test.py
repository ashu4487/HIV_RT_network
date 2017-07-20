# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:39:07 2017

@author: ashutosh
"""

import pickle
with open('all_contacts.txt','rb') as all_con:
    all_contacts=pickle.load(all_con)
with open('frame0contacts.txt','rb') as rcon:
    con_read=pickle.load(rcon)
    
print(len(con_read))
print(len(all_contacts))
#print(con_read)
#print(all_contacts)
#for item in con_read:
#    print(item)
#    if item in list(all_contacts):
#        index=list(all_contacts).index(item)
#        #contact_mat[index,count]=1
#        #count=count+1
#        print(index)