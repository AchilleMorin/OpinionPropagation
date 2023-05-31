# -*- coding: utf-8 -*-
"""
Created on Fri May 19 09:49:32 2023

@author: Achille
"""


{'beta_leeves': 2, 'mu_leeves': 0.1,'mu_inner':1, 'troncate': 6, 'initial fraction': 0.2}





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
        for c in node.children:
            if lattice.lattice[c]==1:
                x+=1/K
            return x*mu_leeves
    
        
    #if node is a leeve
        
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