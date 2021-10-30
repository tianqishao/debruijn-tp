#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
import pickle  



__author__ = "Tianqi Shao"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Tianqi Shao"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Tianqi Shao"
__email__ = "amandineshao135@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()

#######################################
# 1-a:identification des k-mers unique
#######################################
def read_fastq(fastq_file):
    if(isfile(fastq_file)):
        with open(fastq_file,'r') as f:
            for line in f:
                yield next(f)[:-1]
                next(f)
                next(f)
"""
fonction de liser la sequence et la couper en une liste de k-mer
"""
def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]
"""
conscrtuir la dictionaire de kmer
"""
def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for i in read_fastq(fastq_file):
        for j in cut_kmer(i,kmer_size):
            if j not in kmer_dict:
                kmer_dict[j] = 1;
            else:
                kmer_dict[j] += 1;
    return kmer_dict
#######################################
## 1-b: construction de l'arbe de Debruijn
#######################################
"""
fonction de construire une graphe orientee
"""
def build_graph(kmer_dict):
    G = nx.Graph()
    for kmer in kmer_dict.key():
        #eulerian(veritice are (k-1)-mers,edge are k-mers)
        #prefixe: kmer[:-1]
        #suffixe: kmer[1:]
        G.add_edge(kmer[:-1],kmer[1:],weight = kmer_dict[kmer])
    return G
#######################################
## 2.Parcours du graphe de Debruijn
#trois fonctions sont necessaires
#######################################

#2-1:prend en entrée un graphe et retourne une liste de noeuds d’entrée
"""
fonction de obtenir l'ensemble des noeuds d'entree
"""
def get_starting_nodes(G):
    start = []
    for node in G.nodes:
        if not G.in_edges(node):
            start.append(node)
    return start

#2-2:prend en entrée un graphe et retourne une liste de noeuds de sortie
"""
fonction de obtenir l'ensemble des noeuds de sortie
"""
def get_sink_nodes(G):
    end = []
    for node in G.nodes:
        if not G.out_edges(node):
            end.append(node)
    return end
    
#2-3:prend un graphe, une liste de noeuds d’entrée et une liste de sortie et retourne une liste de tuple(contig, longueur du contig)
def get_contigs(G,start,end):
    contigs = []
    for s in start:
        for e in end:
            for path in nx.all_simple_paths(G,s,e):
                cont = path[0]
                for i in enumerate(path):
                    if i == 0:
                        cont = path[i]
                    else: 
                        cont += path[i][-1]
                ##liste de tuple:[cont,len(cont)]
                contigs.append((cont,len(cont)))
    return contigs

def save_contigs(contigs,file_name):
    if(isfile(file_name)):
        with open(file_name,"w") as file:
            for i,cont in enumerate(contigs):
                file.write(">contig_{} len = {}\n".format(i,cont[1]))
                file.write(fill(cont[0]+"\n"))
#######################################
## 3.Simplification du graphe de de bruijn
## a:resolution des bulles
#######################################
"""
Fonction de calculer l'ecart-type
"""
def std(vect_val):
    return statistics.stdev(vect_val)
"""
fonction de calculer le average weight
"""
def path_average_weight(G,path):
    weight = []
    for i in enumerate(path):
        weight.append(G.edges[path[i],path[i+1]]['weight'])
    return statistics.mean(weight)
"""
prend un graphe et une liste de chemin, la variable booléenne
supprimer les paths de graph 
"""

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    ## supprimer les multiples nodes
    for path in path_list:
        # path[entree:sink]
        graph.remove_nodes_from(path[1:-1])
        #supprime le node et les liens adjacents.
        if delete_entry_node:
            graph.remove_node(path[0])#entree
        if delete_sink_node:
            graph.remove(path[-1])#sink
    return graph


"""
fonction de comparer et choisir le meulleur path,puis suppremer les path inutil
"""
def remove_paths(G, path_list, delete_entry_node, delete_sink_node):

	for path in path_list: 
		G.remove_nodes_from(path[1:-1]) 
		if delete_entry_node: 
			G.remove_node(path[0]) 
		if delete_sink_node:
			G.remove_node(path[-1])
	return G	



def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    max_weight = max(weight_avg_list)
    max_length = max(path_length)
    pos_max = list(range(0,len(path_list)))
    del_path = []
    pos_tmp = list(pos_max)
    if std(weight_avg_list)!=0:
        for i in pos_max:
            if weight_avg_list[i]<max_weight:
                del_path.append(path_list)
                pos_max.remove(i)
    path_length = [path_length[i] for i in pos_max]
    pos_tmp = list(pos_max)

    if len(pos_max)>1 and std(path_length)!=0:
        for i in pos_tmp:
            if path_length[i]<max_length:
                del_path.append(path_list[i])
                pos_max.remove(i)
    
    if len(pos_max)>1:
        indice = randint(0,len(pos_max))
        pos_max.remove(indice)
        for i in pos_max:
            del_path.append(path_list[i])
    return remove_paths(graph, del_path, delete_entry_node,delete_sink_node)

"""
fonction détermine les chemins “simples” possible entre deux nœuds (un ancêtre et un descendant) et calculer la longueur et le poids de ces chemins. 
"""
def solve_bubble(graph, ancestor_node, descendant_node):
    bubble_path = []
    bubble_length = []
    bubble_weight = []
    for path in nx.allsimple_path(graph,ancestor_node,descendant_node):
        bubble_path.append(path)
        bubble_length.append(len(path))
        bubble_weight.append(path_average_weight(graph,path))
    return select_best_path(graph,bubble_path,bubble_length,bubble_weight)

"""
prend un graphe “brut” et retourne un graphe sans bulle.
"""
def simplify_bubbles(graph):
    bubble = []
    for node in graph:
        pred_node = list(graph.predecessors(node))
        if len(pred_node)<2:
            anc = nx.lowest_common_ancestor(graph, pred_node[0], pred_node[1]) # le plus petit ancetre commum
            bubble.append([anc,node])
    for i in bubble:
        graph = solve_bubble(graph,i[0],i[1])
    return graph
#######################################
## 3.Simplification du graphe de de bruijn
## b:detection des pointes
#######################################

"""
fonction retire les branches d'entree non optimal
retourne graph sans chemin d'entree indesirables 
"""
def solve_entry_tips(graph, starting_nodes):
    ancestors = []
    path_list = []
    path_weight = []
    path_length = []

    for node in starting_nodes:
        for next in nx.descendants(graph,node):
            if len(graph.pred[next])>=2 and next not in ancestors:
                ancestors.append(next)
    for node in starting_nodes:
        for i in ancestors:
            for path in nx.all_simple_paths(graph,node,i):
                path_list.append(path)
                path_weight.append(path_average_weight(graph,path))
                path_length.append(len(path))
        graph = select_best_path(graph,path_list,path_length,path_weight)

    return graph
"""
supprimer les sink non optimal
"""
def solve_out_tips(graph, ending_nodes):
    descendants = []
    path_list = []
    path_weight = []
    path_length = []
    for node in ending_nodes:
        for next in nx.ancestors(graph,node):
            if len(graph.succ[next])>=2 and next not in descendants:
                descendants.append(next)
    for node in ending_nodes:
        for i in descendants:
            for path in nx.all_simple_paths(graph,node,i):
                path_list.append(path)
                path_weight.append(path_average_weight(graph,path))
                path_length.append(len(path))
        graph = select_best_path(graph,path_list,path_length,path_weight)

    return graph


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    print(args)

    kmer_size = args.kmer_size
    fastq_file = args.fastq_file
    kmer_dict = build_kmer_dict(fastq_file,kmer_size)
    graph = build_graph(kmer_dict)

    start_node = get_starting_nodes(graph)
    sink_node = get_sink_nodes(graph)

    graph = simplify_bubbles(graph)
    
    graph = solve_entry_tips(graph,start_node)
    graph = solve_out_tips(graph,sink_node)

    contig_list = get_contigs(graph,start_node, sink_node)

    save_contigs(contig_list, args.output_file)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
