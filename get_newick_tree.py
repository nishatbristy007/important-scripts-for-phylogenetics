#CREATE newick trees from dots and save in a file 'trees.nwk' for CASet and DISC distances.

import networkx as nx
from Bio import Phylo
import numpy as np
import pandas as pd
from networkx.drawing.nx_agraph import read_dot
from itertools import islice
import os

def to_newick_recursively(G, node, multilabels):
    """Recursively convert a NetworkX tree graph G to Newick format starting at the given node.
    Internal nodes might have multiple labels, which are stored in the multilabels dictionary."""
    
    newick = ""

    
    children = list(G.successors(node))
    if children:
        newick += "("
        for child in children:
            newick += to_newick_recursively(G, child, multilabels) + ","
        newick = newick.rstrip(",")  # Removing trailing comma
        newick += ")"

    # Add the current node's labels
    labels = multilabels.get(node, [])
    if len(labels) ==1:
        label_str = labels[0]
    else:
        label_str = "{" + ",".join(labels) + "}"
    newick += label_str

    return newick

def check_inputs(folder):
    flag=1
    if not os.path.isfile(folder+'/scarlet_output/scarlet_output.B'):
        print('scarlet output not found')
        flag=1
    if not os.path.isfile(folder+'/outputs/medicc2_output/clusters.tsv'):
        print('medicc2 clusters not found')
        flag=0
    
    if not os.path.isfile(folder+'/Z.tsv'):
        print('Z not found')
        flag=0
    if not os.path.isdir(folder+'/outputs/output_sctusvext_2/'):
        print('sctusvext output dir not found')
        flag=0
    if not os.path.isfile(folder+'/outputs/output_sctusvext_2/T.dot'):
        print('sctusvext output dir not found')
        flag=0
    if not os.path.isdir(folder+'/outputs/medicc2_output/'):
        print('medicc2 output dir not found')
        flag=0
    
    return flag

def get_multilabels(Z):
    multilabels = {}
    #print(Z.shape)
    for i in range(Z.shape[1]):
        cells = np.where(Z[:,i]==1)[0]
        
        if len(cells) ==1:
            cells +=1
            cells = ['c'+str(cells[0])]
        elif len(cells)>1:
            cells +=1 # cells are 1 indexed.
            cells = ['c'+str(j) for j in cells] #adding c 
        else:
            cells = ['0']
        multilabels[str(i)] = cells
    display(Z)
    print(multilabels)
    return multilabels

def get_multilabels_est(Z,clusters,true_tree):
    multilabels = {}
    
    for i in range(Z.shape[1]):
        cells = np.where(Z[:,i]==1)[0]
        
        if len(cells) ==1:
            cells +=1
            
            label = clusters[clusters['cluster'] == cells[0]]['cell_no'].iloc[0]
        
            cells = ['c'+str(label)]
        elif len(cells)>1:
            cells +=1 # cells are 1 indexed.
            cells = ['c'+str(j) for j in cells] #adding c 
        else:
            cells = ['0']
        multilabels[str(i)] = cells
    return multilabels

def convert_one_idxd_dot_to_zero_idx(G):
    H = nx.MultiDiGraph()
    for node in G.nodes():
        new_name = str(int(node) - 1)
        H.add_node(new_name)
    # Iterating over the nodes and rename them
    for node in G.nodes():
        for edge in G.out_edges(node):
            H.add_edge(str(int(edge[0]) - 1), str(int(edge[1]) - 1))
        # Remove the original node from the copied graph
    return H



def save_newick(folder, output_file, tree_indexed,clone,true_tree,clusters):
    G = read_dot(folder+'/T.dot')
    Z = np.array(pd.read_csv(folder+'/Z.tsv',sep='\t',header=None))[:clone,:clone]
    if tree_indexed == 1:
        G = convert_one_idxd_dot_to_zero_idx(G)
        multilabels = get_multilabels(Z)
    if true_tree == 0:
        multilabels = get_multilabels_est(Z,clusters,true_tree)
    #print(multilabels)
    root = [node for node, indegree in G.in_degree() if indegree == 0][0]
    newick_str = to_newick_recursively(G, root, multilabels)+';\n'
    
    # Save the Newick string to a file
    print(newick_str)
    with open(folder+'/'+output_file, 'w') as f:
        f.write(newick_str)
    return newick_str

true_tree=0
newick_est = save_newick('sctusvext_output_dir/', 'T_est.nwk',0,clone,true_tree,clusters)
true_tree=1
newick_true = save_newick('Sctusvext_input_parent_dir/', 'T_true.nwk',1,clone,true_tree,clusters)
        
