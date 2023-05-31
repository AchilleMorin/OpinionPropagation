# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 15:10:17 2023

@author: Achille
"""

#%% General imports
import numpy as np
import matplotlib.pyplot as plt

#%% Local Imports



    
#%% DISTANCES AND TOPOLOGY OF THE LATTICE

def test_topo():
    print('test_topo')


def mrca(x,y): #stands for most recent common ancestor
    if x==0 or y==0:
        return 0    
    a=[]
    n=min(len(x),len(y))
    for i in range(n):
        if x[i]==y[i]:
            a.append(x[i])
        else:
            return tuple(a)
    return tuple(a)

def distzero(x):
    if x==0:
        return 0
    return len(x)-1
        

def distance(x,y):
    '''
    x and y have to be tuples or zero
    '''

    if x==0:
        return distzero(y)
    if y==0:
        return distzero(x)
    k=len(x)+len(y)-2*len(mrca(x,y))
    return k


#%% SEARCH IN A TREE


def search_level(lattice,target=1):  #returns the nodes of target length in a list
    '''
    lattice has to be a list; because an object of type Lattice is not iterable.
    '''
    result=[]
    for node in lattice:
        if node==0:
            if target==1:
                result.append(node)
        elif type(node)==tuple:
            if len(node)==target:
                result.append(node)
        else:
            if node.len==target:
                result.append(node)
    return result

def find_leeves(lattice):
    if len(lattice)==1:
        return [0]
    k=0
    for x in lattice:
        if x==0:
            continue
        if len(x)>k:
            k=len(x)
    leeves=search_level(lattice,target=k)
    return leeves

def find_IDs(arbre:dict):
    '''
    arbre doit etre un dictionnaire'''
    zeros=[]
    ones=[]
    twos=[]
    for key in arbre:
        if arbre[key]==0:
            zeros.append(key)
        elif arbre[key]==1:
            ones.append(key)
        elif arbre[key]==2:
            twos.append(key)
    return zeros,ones,twos





#%% FIlIATION in a LATTICE

def descendance(x,lattice):
    
    if x==0:
        return list(lattice.keys())  
    nodes=list(lattice.keys())
    depth=len(nodes[-1])-1             #ameliorer ca
    lenx=distzero(x)
    k=depth-lenx
    des=[x]
    for i in range(k):
        mothers=search_level(des,target=lenx+i+1)
        for mo in mothers:
            chil=children(mo,lattice)
            for c in chil:
                if c not in des:
                    des.append(c)
    return des
        
        
    

def children(x,lattice):
    digit=0
    chil=[]
    X=x
    if X==0:
        while (0,digit) in lattice:
            chil.append((0,digit))
            digit+=1
    else:
        if type(X)!=tuple:
            X=x.coor
            if X==0:
                return chil
        y=X+(digit,)
        while y in lattice:
            chil.append(y)
            digit+=1
            y=X+(digit,)
    return chil


def mother(x):
    
    if x==0:
        #print('Warning: root has no mother')
        return 0
    return x[:-1]

def layer(x,lattice):
    return search_level(lattice,target=distzero(x)+1)
    


def brothers(x,lattice):
    if x==0:
        #print('Warning: root has no brother')
        return 0
    bros=[]
    mo=mother(x)
    lay=layer(x,lattice)
    for i in lay:
        mo_i=mother(i)
        if mo_i==mo:
            bros.append(i)
    return bros
    


