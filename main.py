# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 11:08:15 2017
This is the main execution script for converting trajectories to dynamic networks and doing a bunch of stuff with them.
Files formed- all_contacts.txt, framesXXXX_contacts.txt, comntact_matrix.txt

@author: ashutosh
"""

import os
import sys
import command_args
import create_network
import numpy as np
import matplotlib.pyplot as plt


#Parse the command line arguments
xtc,tpr,ndx,skip,b,e,out=command_args.parseargs()
##### Determining all unique contacts in the each frame 

if os.path.isfile("all_contacts.txt"):
    print("Contacts already existing in all_contacts.txt")
    with open('all_contacts.txt','r') as all_con:
        all_contacts=[tuple(i_ac.strip().split(' ')) for i_ac in all_con]
        #print(all_contacts)

else:
    print("Fresh run")
    all_contacts_set=set()                                      #create empty set to keep all contacts
    num_frames=0
    for i in range(int(b),int(e),int(skip)):
        num_frames= num_frames+1   
        fname=create_network.extract_frame(xtc,tpr,ndx,i) #Extracts frame and saves the PDB file with fname
        contacts=create_network.find_contact(fname) # Find the contacts in protein structure
        outf=open(fname[:-4]+'contacts.txt','a')
        for temp_c in contacts:
            outf.write(' '.join(str(s) for s in temp_c) + '\n')
        outf.close()            
            
        scon=set()
        for con in contacts:
            scon.add(con)
            print(con)
        all_contacts_set=set(all_contacts_set).union(scon) #Add new contacts to all contact list
        #print(list(all_contacts)[0:5])
        #print(con_count)
        print(len(all_contacts_set))
        #print(sys.getsizeof(all_contacts))
        os.remove(fname)                            #removes the PDB file to save space    
    all_contacts=list(all_contacts_set)
    all_out=open('all_contacts.txt','a')
    for temp_ac in all_contacts:
        all_out.write(' '.join(str(s) for s in temp_ac) + '\n')
    all_out.close() 

#### Forming contact matrix
    
if os.path.isfile("contact_matrix.txt"):
    print("Contact matrix already exists in folder. Now reading it")
    contact_mat=np.loadtxt('contact_matrix.txt')
    [num_cont,num_frames]=np.shape(contact_mat)
    print(num_cont,num_frames)
    
else:
#    num_frames=100 #### only for test
    contact_mat=np.zeros([len(all_contacts),num_frames])
    count=0
    for i1 in range(int(b),int(e),int(skip)):
        with open('frame'+str(i1)+'contacts.txt','r') as rcon:
            con_read=[tuple(i_c.strip().split(' ')) for i_c in rcon]
            print(con_read)
        for item in con_read:
            if item in all_contacts:
                index=all_contacts.index(item)
                contact_mat[index,count]=1
        count=count+1
    np.savetxt('contact_matrix.txt',contact_mat,fmt="%d")
    
#plt.imshow(contact_mat,aspect="auto",interpolation='nearest')
#[num_cont,num_frames]=np.shape(contact_mat)

sum_contact=np.sum(contact_mat,axis=1)
cont_prob=sum_contact/num_frames    #Calculate contact probability
print(sum_contact)
print(cont_prob)
# print(np.shape(cont_prob))
 
contacts_selected,sel_cont_rows=create_network.extract_dynamic_contacts(0.2,0.7, all_contacts,contact_mat,num_frames)  
print(len(contacts_selected))
#print(len(sel_cont_rows))
print(np.shape(sel_cont_rows))
create_network.PCA_contact_mat(sel_cont_rows)
#plt.plot(range(len(cont_prob)),sorted(cont_prob,reverse=True))
#plt.show() 