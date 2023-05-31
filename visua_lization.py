# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:22:53 2023

@author: Achille
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 16:12:29 2023

@author: Achille
"""
#%% General imports
import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
import glob
from PIL import Image
import os


#%% Local Imports

os.chdir('C:/Users/Achille/Documents/ENS/M2/Stage/Code/arbre_fini/3_etats')


from topo_logy import find_leeves,find_IDs,search_level


from some_tools import conn_dic,hierarchy_pos,order,length




#%% Visualizing the leeves of a lattice

        
def plot_leeves(lattice,title=None,save=False,show=True):
    
                            #il faut trouver mieux ici. il se pourrait que nos noeuds soient desordonnes. Donc regarder le max
    leeves=find_leeves(lattice)       #mais on ne peut pas regarder len(0) (noeud origine)
    n=len(leeves)                           #faire de la POO pour definir proprement ce que c'est qu'un noeud, qu'un arbre...
    plt.figure(figsize=(n,3))
    ax1 = plt.axes()
    ax1.axes.get_yaxis().set_visible(False)
    for i in range(n):
        if lattice[leeves[i]]==1:
            plt.scatter(i,1,color='red',marker='*')
        else:
            plt.scatter(i,0,color='blue',marker='.')
    if title:
        plt.title(title)
        if save:
            plt.savefig(title)
    if show:
        plt.show()
        
#%% CIRCULAR PLOTS OF A WHOLE LATTICE


def circular_plot(lattice,highlight=[],highlight2=[],root=0,save=False,filename=None,title=None,circular=True):
    
    CD=conn_dic(lattice)
    G=nx.Graph(CD)
    pos=hierarchy_pos(G,root=root,width = 2*np.pi, xcenter=0)
    if circular:
        new_pos = {u:(r*np.cos(theta),r*np.sin(theta)) for u, (theta, r) in pos.items()}
    else:
        new_pos=pos
    plt.figure()
    nx.draw(G,pos=new_pos,node_size=5)
    nx.draw_networkx_nodes(G, pos=new_pos, nodelist = highlight, node_color = 'red', node_size = 50)
    if highlight2:
        nx.draw_networkx_nodes(G, pos=new_pos, nodelist = highlight2, node_color = 'green', node_size = 50)
    plt.title(title)
    if save:
        plt.savefig(filename,bbox_inches='tight')
    pass

#%% BUILDING VISUALIZATION FOR CONTINUOUS TIME TRAJECTORIES


def continuous_plot(states,jumps,nb_frames=400,circular=True):
    
    root=0                    #do better: root(lattice),lattice.root
    jumps=np.array(jumps)
    frame=states[0]
    CD=conn_dic(frame)   #resumer tout ça en une fonction
    G=nx.Graph(CD)
    pos=hierarchy_pos(G,root=root,width = 2*np.pi, xcenter=0)
    if circular:
        new_pos = {u:(r*np.cos(theta),r*np.sin(theta)) for u, (theta, r) in pos.items()}
    else:
        new_pos=pos
    
    tmax=jumps[-1]
    dt=tmax/nb_frames
    
        
    for k in range(nb_frames+1):
        
        i=np.max((np.where(jumps<=k*dt)))
        identities=find_IDs(states[i])
        highlight,highlight2=identities[1],identities[2]
        plt.figure()
        nx.draw(G,pos=new_pos,node_size=5)
        nx.draw_networkx_nodes(G, pos=new_pos, nodelist = highlight, node_color = 'red', node_size = 50)
        nx.draw_networkx_nodes(G, pos=new_pos, nodelist = highlight2, node_color = 'green', node_size = 50)
        plt.title(f't={k*dt}')
        plt.savefig("{0:05}".format(k),bbox_inches='tight')
        plt.close()
    return
        

#%% VISUALIZING THE FRACTION OF INFECTED THROUGH TIME

def time_follow_x_layers(states,jumps,lattice):
    
    K=lattice.order
    tmax=jumps[-1]
    jumps=np.array(jumps)
    T=len(jumps)
    
    if lattice.root==0:
        n=lattice.len
    else:
        n=lattice.len+1
    
    data=np.zeros([T,n])
      
        
    for t in range(T):
        
        for index in range(n):
            
            layer=getattr(lattice, f'layer{index}')
            K=len(layer)
            x=0
            for node in layer:
                if states[t][node]==1:
                    x+=1
            x*=1/K 
            
            data[t,index]=x
        
    
    #totfrac=sum(data[:,i]*len(layers[i]) for i in layers)
    #totfrac*=1/N
    jumps=list(jumps)
    jumps.append(jumps[-1])
    plt.figure()
    for i in range(n):
        if i==0:
            plt.scatter(jumps[:-1],data[:,i],label=f'Layer {i}')
        else:
            plt.stairs(data[:,i],jumps,label=f'Layer {i}')
    #plt.plot(jumps,totfrac,label='Total fraction',color='black',linestyle='dotted')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Convinced fraction pro layer')
    plt.title('Convinced fractions in each layer')
    plt.show()

    return 


#%%



def time_follow_x_fratries(states:list,jumps:list,lattice,nb_frames=400):
    
    aunts=lattice.root
    K=lattice.order
    tmax=jumps[-1]
    jumps=np.array(jumps)
    dt=tmax/nb_frames
    data=np.zeros([nb_frames,1+len(aunts)])
    times=[]
        
    for t in range(nb_frames):
        times.append(t)
        index=np.max(np.where(jumps<=t*dt))
        data[t,-1]=jumps[index]
        for m,aunt in enumerate(aunts):
            
            node_aunt=lattice.nodes[aunt]
            fraction=0
            for node in lattice.leeves:
                if node in node_aunt.descent:
                    fraction+=states[index][node]
            fraction*=1/K
            data[t,m]=fraction
    
    plt.figure()
    for m in range(len(aunts)):
        plt.plot(data[:,-1],data[:,m],label=f'Family {m}')
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Convinced fraction pro Familly')
        plt.title('Convinced fractions in each family')
    plt.show()
    
        
    
    
    
    
    
    
    


        
 
#%% GROUPING TREES INTO A GIF

def make_gif(path=None,duration=100):
    if path:
        os.chdir(path)
    frames = []
    imgs = glob.glob("*.png")
    list.sort(imgs)                    #pour l'instant ca marche mais il faut ordonner les images (2 vient entre 19 et 20)
    for i in imgs:
        new_frame = Image.open(i)    
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save('png_to_gif.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=duration) #loop=0