#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:48:17 2020

@author: nehabinish
"""

import numpy as np
import time


start = time.time()

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
            # creating lists for the linked lists
            i.append(int(c[0])) # out node
            j.append(int(c[1])) # in node

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
                
    nodes = sorted(nodes) # sorted values of nodes
    
    N = len(nodes) # total number of nodes        
     
    return nodes,N
        

# finding the linked nodes to define adjacent matrix  
link_i,link_j = read_file('adjacency.txt')

# calculating no of nodes and getting a list of nodes
nodes,N = Nodes(link_i,link_j)

# DIRECTED NETWORK
A = np.zeros((N,N))

# adjacent matrix
for i in range(0,len(link_i)):
    A[link_j[i]-1][link_i[i]-1] = 1


print('')
print('Adjacency Matrix')               
print(A)

# kj_out - out degree of the jth node
kj_out = np.zeros(N)
for i in range(0,N):
    kj_out[i] = sum(A[:,i])

print('')
print('kj_out')
print(kj_out)

# Stochastic matrix
S = np.zeros((N,N))

for i in range(0,N):
    for j in range (0,N):
        
        if (kj_out[i]==0):
            S[j][i] = 1/N
        else:
            S[j][i] = A[j][i]/kj_out[i]

print('')
print('Stochastic matrix') 
print(S)   
      
alpha = 0.85

# first preferential vector v with every element having equal probabilities
v0 = np.zeros(N)
for i in range(0,N):
    v0[:] = 1/N
        
# Perron-Frobenius operator G  
G = np.zeros((N,N))
for j in range(0,N):
    for i in range(0,N):
        
        G[i][j] = S[i][j]*alpha + (1-alpha)*v0[i]

print('')
print('Perron-Frobenius operator') 
print(G) 

# finding the final probability vector
P = []
P.append(v0)

v1 = np.dot(G,v0)

# defining episilon
epsilon = np.zeros(N)
for i in range(0,N):
    epsilon[:] = 0.000001
 
# finding the steady state probability distribution vector P
for i in range(0,1000):
    P.append(np.dot(G,P[i]))

    if(np.all((P[i+1]-P[i]) < epsilon)):
        pagerank_vector = P[i+1]
        print('')
        print('no of iterations: {}'.format(i))
        break;

print('')
print('Page-rank Vector') 
print(pagerank_vector)

# checking if the page rank vector is correct
print('')
print("Total probability of the page rank vector")
print(sum(pagerank_vector))

# to make a dictionary with nodes and probabilities
# using zip()
print('') 
print('Nodes and their probabilities')
pg_rnk = dict(zip(nodes,pagerank_vector)) 
print(pg_rnk)

# to show ranking of nodes
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
    
    
# k_in - in degree of the node
k_in = np.zeros(N)
for i in range(0,N):
    k_in[i] = sum(A[i,:])

    
print('')
print('k_in')
print(k_in)

print('It took {0:0.1f} seconds'.format(time.time() - start))