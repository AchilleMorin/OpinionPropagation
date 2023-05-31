# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:04:58 2023

@author: Achille
"""
#%% General imports
import numpy as np
import os

#%% Local Imports

os.chdir('C:/Users/Achille/Documents/ENS/M2/Stage/Code/arbre_fini/3_etats')
from para_meters import parameters
from topo_logy import search_level  

#%% CONSTRUCTION OF AN EMPTY LATTICE


def essai():
    return 'OUI'
    
    
    
def daughters(node,order=2): #returns the daughters of a given node
    dghtrs=[]
    if node==0:
        return [(node,i) for i in range(order)]
    for i in range(order):
        dghtrs.append(node+(i,))
    return dghtrs


def make_ltc(n,order=2,root=[0]):  
    ''' Creates a lattice. k is the distance between the root and the leeves.
    Warning: if the root is other than 0 (for example for two-families trees
    where the root should be given as [(0,0),(0,1)], n will be the distance
    between (0,0) and the leeves. So the total length of the tree is n+1.'''
    
    lattice=[0]
    if root==[0]:
        size=1
    else:
        size=len(root[0])
    for i in range(n):
        for node in root:
            if node!=0:
                lattice.append(node)
            newnodes=daughters(node,order=order)
            for nn in newnodes:
                lattice.append(nn)
        root=search_level(lattice,target=i+size+1)   #i+2
    lattice={node:0 for node in lattice}
    return lattice


#%% INITIAL DENSITIES:
'''ces fonctions proposent différentes façons d'initialiser la contamination
    dans une population'''
    
def single_contaminated(n,params=None):
    '''une seule est initialement contaminée'''
    out=[0 for i in range(n)]
    out[0]=1
    return out

def initial_fraction(n,params=parameters):
    '''une fraction donnée de feuilles aléatoires est contaminée'''
    fraction=params['initial fraction']
    out=[0 for i in range(n)]
    number=int(fraction*n)
    ones=np.random.choice(np.arange(n),size=number,replace=False)
    for i in ones:
        out[i]=1
    return out

def single_contaminated_aunts(n,aunts=2,params=None):
    return
    
    
    

