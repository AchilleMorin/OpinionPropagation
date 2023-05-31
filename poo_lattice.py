# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:58:32 2023



Ces arbres sont construits avec des individus pouvant avoir trois
identités diférentes:
    les 0 ne sont pas contaminés, ni contaminants
    les 1 sont contaminés et contaminants
    les 2 ont été contaminés par leur mère, mais ils ne sont pas contaminants

@author: Achille
"""

#%% IMPORTS

import os
import time
import numpy as np

#%% LOCAL IMPORTS

os.chdir('C:/Users/Achille/Documents/ENS/M2/Stage/Code/arbre_fini/3_etats')

from visua_lization import circular_plot,continuous_plot,time_follow_x_fratries,time_follow_x_layers,make_gif
from para_meters import parameters
from conta_mination import trajectory  #Warning: circular imports
from topo_logy import find_leeves,distzero,mother,brothers,children,descendance
from topo_logy import search_level,test_topo
from build_lattice import single_contaminated,initial_fraction,make_ltc,essai
from some_tools import conn_dic, estula


#%% DEFINITION OF A LATTICE

class Lattice:
    '''
    This class creates a lattice object that has many attributes.
    
    For a regular lattice, just enter n: total length of the lattice and K: order of the lattice.
    
    Root can be changed manually from [O] (default) to [(0,0),(0,1)] for example for two-families lattices.
    
    More generally, for a tree with m aunts and m families, set aunts=m.
    Otherwise, aunts should always be set to 0 (default).
    '''
    #les attributs devraient etre des choses qui ne changent jamais
    def __init__(self,n:int,K:int,root=[0],aunts=0,display_progression=True):  #ajouter ROOT!!!!!!!!!!!!!!
        
        t0=time.time()
        
        if aunts!=0:
            root=[(0,i) for i in range(aunts)]
        
        #strucutre of the lattice
        self.lattice=make_ltc(n,K,root=root)
        self.order=K
        self.len=n
        self.root=root
        self.size=len(self.lattice)
        
        if root==[0]:
            self.type='regular'
        else:
            self.type='several roots. Warning: the total length of the tree is n+1; n being the distance between the non-0 roots and the leeves.'
        
        
        #LAYERS
        self.leeves=find_leeves(self.lattice) 
        
        progress=[]
        
        n_layers=self.len+1
        if aunts!=0:
            n_layers+=1
        for layer in range(n_layers):
            setattr(self, f'layer{layer}', 
                    search_level(self.lattice,target=layer+1)) #en faire un dico avec les id?
            
            
            if display_progression:
                perc=int(100*layer/(self.len+1))
                if perc//10 not in progress:
                    progress.append(perc)
                    turn=int((time.time()-t0)//60)
                    print(f'Building layers. Completed: {perc}%. Runtime: {turn} minutes')
            
        
        #NODES
        progress=[]
        
        self.nodes={}
        for (i,node) in enumerate(self.lattice):
            self.nodes[node]=Node(node,self.lattice)
            
            if display_progression:
                perc=int(100*i/self.size)
                if perc//1 not in progress:
                    progress.append(perc)
                    turn=int((time.time()-t0)//60)
                    print(f'Building nodes. Completed: {perc}%. Runtime: {turn} minutes')
            
            
    def codic(self):
        return conn_dic(self.lattice)
        
        
    def zeros(self):
        z=[]
        for node in self.nodes:
            if self.nodes[node].id == 0:
                z.append(node)
        return z
    
    def ones(self):
        o=[]
        for node in self.nodes:
            if self.nodes[node].id==1:
                o.append(node)
        return o
    
    def twos(self):
        t=[]
        for node in self.nodes:
            if self.nodes[node].id==2:
                t.append(node)
        return t
    
    
    def viz(self,circular=True,highlight='root',highlight2=None):
        if highlight=='root':
            circular_plot(self.lattice,circular=circular,highlight=self.root)
        elif highlight=='IDs':
            circular_plot(self.lattice,circular=circular,highlight=self.ones(),highlight2=self.twos())        
        elif highlight==None:
            circular_plot(self.lattice,circular=circular)
        elif type(highlight)==list:
            circular_plot(self.lattice,circular=circular,highlight=highlight,highlight2=highlight2)
        

            
    
    def __repr__(self):
        return f'Lattice with length {self.len}, order {self.order} and root {self.root}'
    
        
    def initial_ID_leeves(self,method=single_contaminated,params=parameters,
                          each_family=False):
        '''This method starts a contamination from the leeves, by setting
        some leeves with ID 1, accroding to the given method. See methods for
        help:
            single_contaminated
            initial_fraction'''
        
        for i,aunt in enumerate(self.root): #root is 0 or [(0,i) for i in range(aunts)]
    
            if not each_family:
                if i>0:
                    break
            
            node_aunt=self.nodes[aunt]
            des=node_aunt.descent  
            aunt_leeves=find_leeves(des)
            distr=method(len(aunt_leeves),params=parameters) 
    
            for i,node in enumerate(aunt_leeves):
                self.nodes[node].id=distr[i]
                self.lattice[node]=distr[i]
                
        
    def clean_ID(self):
        for node in self.nodes:
            self.nodes[node].id=0
            self.lattice[node]=0
            
        
        
        
            
            
            
#%% DEFINITION OF A NODE            

class Node:
    
    def __init__(self,coordinates:tuple,lattice):
        self.coor=coordinates
        self.distzero=distzero(self.coor)
        self.mother=mother(self.coor)                     #est-ce que tout ça ça doit être des coordonnées
        self.brothers=brothers(self.coor,lattice)        #ou des vrais noeuds
        self.children=children(self.coor,lattice)     #si c'est des noeuds on va jamais s'en sortir nan, trop récursif
        self.descent=descendance(self.coor,lattice)
        self.id=0
        self.cousins=search_level(lattice,target=self.distzero+1)
        if self.coor==0:
            self.len=1
        else:
            self.len=len(self.coor)
        
        
        
        
        def layer(x,lattice):
            return search_level(lattice,target=distzero(x)+1)
        
        
    def __repr__(self):
        return f'Node with coordinates {self.coor} and identity {self.id}'
    


