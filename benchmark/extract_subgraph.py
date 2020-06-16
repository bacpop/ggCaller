"""
Python scripts that, froma huge graph file, extract a sequence of smaller subgraph
"""
import numpy as np
import queue

from networkx.algorithms.traversal.breadth_first_search import bfs_edges

def random_source(gfa_):
    nodes = gfa_.nodes()
    casual_index = np.random.randint(0, len(nodes)-1)

    return nodes[casual_index]

def bfs_custom(gfa_, source, distance_source):
    nodes_dict = {}
    node_queue = queue.Queue()
    node_queue.put(source)
    nodes_dict[source] = 1

    while not node_queue.empty():
        node = node_queue.get()

        node_neighbors = gfa_.neighbors(node)
        for neighbor in node_neighbors:
            if not nodes_dict.get(neighbor):
                nodes_dict[neighbor] = nodes_dict.get(node) + 1
                if nodes_dict.get(neighbor) - nodes_dict.get(source) <= distance_source:
                    node_queue.put(neighbor)
                else:
                    nodes_dict[neighbor] = False

    return nodes_dict

def print_file_nodes(nodes_dict, source, number_source, grade, distance):
    
    number_source = str(number_source)
    while not len(number_source) == 3:
        number_source = "0"+number_source
    
    path = 'benchmark/benchmark_graphs/source{}_g{}.txt'.format(number_source, grade)
    with open(path, 'w') as nodes_file:
        i = 0
        if grade == 0:
            nodes_file.write(source+'|')
            i += 1  
        for node in nodes_dict:
            if nodes_dict.get(node) > distance/2 + 1:
                i += 1
                nodes_file.write(node+'|')
                if i >= 500:
                    nodes_file.write('\n ')
                    i = 0
        
        nodes_file.write('\n')
        nodes_file.close()

    return path

def extract_subgraph(gfa_, number_source, max_distance):
    
    with open('benchmark/benchmark_graphs/list_file_node.txt', 'w') as list_file:
        for i in range(0, number_source):
            source = random_source(gfa_)
            for grade in range(0, max_distance+2):
                if grade == max_distance + 1:
                    n_source = str(i)
                    while not len(n_source) == 3:
                        n_source = "0"+n_source
                    path = 'benchmark/benchmark_graphs/source{}_g{}.txt'.format(n_source, grade)
                    with open(path, 'w') as nodes_file:
                        nodes_file.write('\n')
                        nodes_file.close()
                    list_file.write(path.lstrip('benchmark').lstrip('/')+'\n')
                else:
                    distance = 2**grade
                    nodes_dict = bfs_custom(gfa_, source, distance)
                    path = print_file_nodes(nodes_dict, source, i, grade, distance)
                    list_file.write(path.lstrip('benchmark').lstrip('/')+'\n')
        list_file.close()