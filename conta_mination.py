# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:15:51 2023


THIS SCRIPT WORKS WITH THE POO FILE poolattice.py



@author: Achille
"""





#%% General imports

import numpy as np
import os
import time

#%% Local Imports
os.chdir('C:/Users/Achille/Documents/ENS/M2/Stage/Code/arbre_fini/3_etats')

from para_meters import parameters
from topo_logy import distance,layer,find_IDs
from some_tools import choice
#from poolattice2 import Lattice,Node #warning: circular import

#%% COMPUTING CHANGE RATES

def qij(node,lattice,ones:list,params=parameters):
    '''node has to be of Node class
    lattice has to be of Lattice class
    '''

    
    beta,mu,troncate=params['beta'],params['mu'],params['troncate']
    
    K=lattice.order
    chil=node.children
    #x=sum([lattice.lattice[c] for c in chil])/K
    x=0
    for c in chil:
        if lattice.lattice[c]==1:
            x+=1/K
      
    # computing contaminations by the cousins
    if node.coor==0: #lattice.root has no cousins
        return x*mu
    q=x*mu
    cousins=layer(node.coor,ones)
    
    for cous in cousins:
        if type(cous)!=tuple:
            cous=cous.coor
        if cous==node.coor:
            continue    #on ne se contamine pas soi-meme
        d=distance(cous,node.coor)//2
        if d>troncate:                  #règle modifiable selon la distance
            continue
        q+=1/(beta**d) 
        
    return q


#%%
def qij_layer_based(node,lattice,ones:list,params=parameters):
    '''node has to be of Node class
    lattice has to be of Lattice class
    '''
    
    beta_leeves,mu_leeves,mu_inner,troncate=params['beta_leeves'],params['mu_leeves'],params['mu_inner'],params['troncate']
    K=lattice.order
    len_tree=lattice.len + len(lattice.root)-1
    
    
    
    # if node is in the inner tree
    if node.len < len_tree:
        x=0
        for c in node.children:
            if lattice.lattice[c]==1:
                x+=1/K
        return x*mu_inner
        
    
    # if node is a mother
    if node.len == len_tree:
        x=0
        for c in node.children:
            if lattice.lattice[c]==1:
                x+=1/K

        return x*mu_leeves
        
    
        
    #if node is a leeve
    q=0
        
    if node.len == len_tree +1:
        cousins=layer(node.coor,ones)
        for cous in cousins:
            if cous==node.coor:
                continue    #on ne se contamine pas soi-meme
            d=distance(cous,node.coor)//2
            if d>troncate:                  #règle modifiable selon la distance
                continue
            q+=1/(beta_leeves**d) 
        return q

#%% SIMULATING A TRAJECTORY


def trajectory(lattice,params=parameters,tmax=50,stop_mother=True,
               display_progression=True,layer_based=False):
    
    '''
    lattice has to be of Lattice class
    returns a simulation of a markov time continuous process in the form
    of discrete states (Xjn)n and the time of exponential jumps.
    L'option stop_mother fait s'arrêter la trajectoire seulement quand la 
    mère est contaminée'''
    
    #initiation
    jumps=[0]
    l=lattice.lattice.copy()
    states=[l]
    t=0
    #runtime monitoring
    progress=[]
    t_init=time.time()
    
    if lattice.ones()==[]:
        print('Attention: aucun individu contaminé initialement')
        return None,None
    
    
    
    ## vérifier que l'état initial ne nécessite pas de contaminations 
    ## pour la descendance; sinon remplacer l'état initial
    
    
    if stop_mother:
        tmax=np.infty
    
    
    while t<tmax:
        
        #progression
        zeros,ones,twos=find_IDs(lattice.lattice)
        
        if display_progression:
            perc=int(100*len(ones)/lattice.size)
            if perc%1==0:
                if perc not in progress:
                    progress.append(perc)
                    turn=int((time.time()-t_init)//60)
                    print(f'Completed: {perc}%. Runtime: {turn} minutes')
            
        
        if zeros==[]: #contamination is over
            break
        
        liste_qij=[]
        for ind in zeros: #ind is a tupple
            node=lattice.nodes[ind] #node is a Node
            
            # CHOIX DE LA REGLE DE CONTAMINATION
            if layer_based:
                liste_qij.append(qij_layer_based(node,lattice,ones,params=params))
            else:
                liste_qij.append(qij(node,lattice,ones,params=params))
            
        qi=sum(liste_qij)
        if qi==0:
            if display_progression:
                print('Contamination no longer possible. Check maybe parameters beta and mu.')
            return states,jumps
        liste_pij=[qij/qi for qij in liste_qij]
        i=choice(liste_pij)
        
        infected=lattice.nodes[zeros[i]]  #infected is a Node
        
        #contaminating the descendance
        des=infected.descent
        lattice.lattice[infected.coor]=1
        lattice.nodes[infected.coor].id=1
        for coor in des: #coor is a tupple
            if lattice.lattice[coor]==0:
                lattice.lattice[coor]=2
                lattice.nodes[coor].id=2
            
        
        l=lattice.lattice.copy()
        states.append(l)
        temps=np.random.exponential(1/(qi/len(infected.cousins)))  #essayer 1/(qi/K) . l'original est #np.exponential(1/qi)
        t+=temps
        jumps.append(t)
        
    return states,jumps  

