# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:04:58 2023

@author: Achille
"""
#%% GENERAL IMPORTS

import numpy as np
import networkx as nx
import random
import time
import os

#%% LOCAL IMPORTS
os.chdir('C:/Users/Achille/Documents/ENS/M2/Stage/Code/arbre_fini/3_etats')


from topo_logy import children


#%% CARACTERISTIC PARAMETERS OF A TREE


def estula():
    print('YES')

def order(lattice):
    if len(lattice)==1:
        print('Warning: only root in lattice')
        return 0
    k=len(children(0,lattice))
    return k

def length(lattice):  #faire mieux que ça dans la POO: n et K sont données donc on les stocke direct
    
    if list(lattice.keys())==[0]:
        return 0
    max=0
    for i in lattice:
        if i==0:
            continue
        if len(i)>max:
            max=len(i)
    return max-1

#%% CONVERTING A LATTICE INTO A CONNECTION DICTIONNARY


def conn_dic(lattice):
    dico={}
    for x in lattice:
        chil=children(x,lattice)
        if chil:
            dico[x]=chil
    return dico

#%% SIMULATING A RANDOM CHOICE GIVEN A LIST OF PROBABILITIES

def choice(liste_pi):
    cumulated=np.cumsum(liste_pi)
    draw=np.random.random()
    i=0
    while draw>cumulated[i]:
        i+=1
    return i  

#%% COMPUTING PLOTTING POSITIONS FOR A HIERARCHICAL LATTICE

def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5):
    '''this function computes the positions of the nodes in order to have a nice circular plot of a given hierarchical graph
    the root should be precised as root=0'''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')
    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))
            
    def _hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):

        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos

    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)

