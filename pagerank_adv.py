#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 22:55:56 2021

@author: nehabinish
"""

'''
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from time import process_time
from collections import defaultdict
'''

import numpy as np
import numba
from numba import jit
import numpy as np
import scipy as sp
import scipy.linalg as LA



# reading both column in the file
def read_file(filename):
    '''
    --------------------------------------------------------------------------
    function to read the file and store data in the file into two lists 
    giving the linked nodes.
    -------------------------------------------------------------------------- 
    '''
    with open(filename,'r') as data:
        i = []
        j = []
        
        for line in data:
            c = line.rstrip( '\n' ) 
            c = line.split()
            #creating lists for the linked lists
            i.append(int(c[0])) #out node
            j.append(int(c[1])) #in node

    return i,j


def Nodes(C1,C2):
    '''
    --------------------------------------------------------------------------
    function that checks for all the unique values in the list to give
    total number of nodes.
    --------------------------------------------------------------------------    
    '''
    nodes = []
        
    for i in C1:               
        if(i not in nodes):
            nodes.append(i)
    
    for j in C2: 
        if(j not in nodes):
            nodes.append(j)
                
    nodes = sorted(nodes) #sorted values of nodes
    
    N = len(nodes) #total number of nodes        
     
    return nodes,N
        

# Finding the linked nodes to define adjacent matrix  
link_i,link_j = read_file('adjacency.txt')

#Calculating no of nodes and getting a list of nodes
nodes,N = Nodes(link_i,link_j)

k_out = []
k_in = []

for i in nodes:
    k_out.append(link_i.count(i))

'''
damping factor
'''

alpha = 0.85


'''
Creating adjacency list using dictionaries

'''
                                  
# Creating a dictionary containing nodes and their incoming links:
# fill the dict with empty lists for each node
in_dict =  {node: [] for node in nodes}
for i in nodes:
    for j in range(0,len(link_j)):
        if (link_j[j] == i):
            # append link to the list for the node
            in_dict[i].append(link_i[j])
                                                                 
           
adj_list = []

for i in range (0,N):
    if i in in_dict:
       adj_list.append(in_dict[i])
    
                
           
# initialising page rank of the vectors
pg_rnk0 = []
for i in nodes:
    pg_rnk0.append(1/N)

pg_rnk1 = []
for i in nodes:
    pg_rnk1.append(0)
    
#Finding the final probability vector
pgR = []
pgR.append(pg_rnk0) 


#defining episilon
epsilon = []
max_iter = 10000

for i in range(0,N):
    epsilon.append(10**(-8))


for i in range(0,max_iter): 
    sub=[]
    
    for j in in_dict:
     
        sum_ = 0
        for k in in_dict[j]:
                    
            sum_ += (pgR[i][k-1]/k_out[k-1])
        
         # checking for dangling nodes   
        for m in range(0,N):
            
            if k_out[m] == 0: 
                sum_ += (pgR[i][m]*(1/N))
            
            
        sub_prev = (1-alpha)/N + alpha*sum_        
        sub.append(sub_prev)
    
    pgr_prev=sub
    for n in range(0,N):
        pgr_prev[n]= sub[n]/LA.norm(sub,1)  
        
    pgR.append(pgr_prev)        
    
    
    difference = []
    zip_object_diff = zip(pgR[i+1], pgR[i])
    for a,b in zip_object_diff:
        difference.append(abs(a-b))
    
    if(np.all(difference < epsilon)):
        pagerank_vector = pgR[i+1]
        print('')
        print('no of iterations: {}'.format(i))
        break; 
    
     
pagerank_vector = pgR[i+1] 
       
print('')
print('Page-rank Vector') 
print(pagerank_vector)

#Checking if the page rank vector is correct
print('')
print("Total probability of the page rank vector")
print(sum(pagerank_vector))

# to make a dictionary with nodes and probabilities
# using zip()
print('') 
print('Nodes and their probabilities')
pg_rnk = dict(zip(nodes,pagerank_vector)) 
print(pg_rnk)      


#to show ranking of nodes
print('') 
print('The order of ranking of the nodes with their given probabilities')
'''
-------------------------------------------------------------------------------
sorting the values of the dictionary in descending order.
dict.items() - return a tuple of all the items in the dictionary
key=lambda x:x[1] - function that sorts the values
reverse=True - default order is ascending, setting reverse=TRUE for descending
-------------------------------------------------------------------------------
'''
sorted_pgrnk = sorted(pg_rnk.items(), key=lambda x:x[1], reverse=True) 
for i in sorted_pgrnk:
    print(i[0],i[1])