#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 13:55:42 2021

@author: nehabinish
"""

import numpy as np
import scipy.linalg as LA
import time

'''
starting counter for time

'''

start = time.time()


def read_file(filename):
    '''
    --------------------------------------------------------------------------
    function to read the file and store data in the file into two lists
    giving the linked nodes.
    --------------------------------------------------------------------------
    '''
    with open(filename,'r') as data: # open file
        i = [] # create empty list to hold column 1
        j = [] # create empty list to hold column 2

        for line in data:
            c = line.rstrip( '\n' ) # removing any unwanted spaces
            c = line.split()        # spliting the text file into two columns

            #creating lists for the linked lists
            i.append(int(c[0]))     # outgoing links
            j.append(int(c[1]))     # incoming links

    return i,j


# Finding tthe outgoing and incoming links
link_i,link_j = read_file('dawiki.txt')

N = max(link_j) # Total number of nodes

alpha = 0.85 # Damping factor

# Creating a range for nodes
nodes = []
for i in range(0,N):
    nodes.append(i+1)

'''
Total number of outgoing links for each node

'''
k_out = []

for i in nodes:
    k_out.append(link_i.count(i))


'''
Creating adjacency list

'''
adj_list = []
for i in nodes:
    inlink = []
    for j in range(0,len(link_j)):
        if (link_j[j]==i):
            inlink.append(link_i[j])

    adj_list.append(inlink)


# Initialising page rank of the vectors
pg_rnk0 = []
for i in nodes:
    pg_rnk0.append(1/N) # Every node has equal probability


# Creating a list to store all page rank vectors
pgR = []
pgR.append(pg_rnk0)


'''
We find the page rank vector either after 100 iterations or with a tolerance of epsilon

'''
# Defining episilon - breaking condition
epsilon = []
max_iter = 100 # Maximum number of iterations

for i in range(0,N):
    epsilon.append(10**(-8))


for i in range(0,max_iter):
    sub=[]

    for j in range(0,N):

        sum_ = 0

        # Checking all the incoming links and augumenting the sum
        for k in adj_list[j]:

            sum_ += pgR[i][k-1]/k_out[k-1]

         # Checking for dangling nodes and augumenting the sum
        for m in range(0,N):

            if k_out[m] == 0:
                sum_ += pgR[i][m]*(1/N)

        # Augumenting by the damping factor
        sub_prev = (1-alpha)/N + alpha*sum_
        sub.append(sub_prev)

    pgr_prev=sub

    for n in range(0,N):
        pgr_prev[n]= sub[n]/LA.norm(sub,1) # Gp(k-1)/||Gp(k-1)||

    pgR.append(pgr_prev) # Storing the previous page ranks


    # Checking if the page ranks follow a tolerence of epsilon
    difference = []
    zip_object_diff = zip(pgR[i+1], pgR[i])
    for a,b in zip_object_diff:
        difference.append(abs(a-b))

    if(np.all(difference < epsilon)):
        pagerank_vector = pgR[i+1]
        print('')
        print('no of iterations: {}'.format(i))
        break; 


pagerank_vector = pgR[i+1] # Final page rank vector

Total_prob = sum(pagerank_vector) # Finding total probability of the page rank vectors

# To make a dictionary with nodes and probabilities
pg_rnk = dict(zip(nodes,pagerank_vector))

'''
---------------------------------------------------------------------------
sorting the values of the dictionary in descending order.
dict.items() - return a tuple of all the items in the dictionary
key=lambda x:x[1] - function that sorts the values
reverse=True - default order is ascending, setting reverse=TRUE for descending
-------------------------------------------------------------------------------
'''
sorted_pgrnk = sorted(pg_rnk.items(), key=lambda x:x[1], reverse=True)

# Writing the sorted page vector into text file
with open('dawiki_sol.txt', 'w') as f:
    print(sorted_pgrnk, file=f)

# Computation time
print('Computation Time - {0:0.1f} seconds'.format(time.time() - start))
