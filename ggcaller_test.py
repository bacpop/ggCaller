#imports
from pygfa import *
import networkx as nx
import re
from Bio.Seq import Seq

#generate graph
#test_graph = pygfa.gfa.GFA.from_file("test.gfa")
#test_graph_3 = pygfa.gfa.GFA.from_file("test3.gfa")
#group3_graph = pygfa.gfa.GFA.from_file("group3_SP_capsular_gene_bifrost.gfa")

#graph = group3_graph

#length graph
#len_graph = len(graph._graph.node)

#takes tsv file from Bifrost query
def add_colours(colours_tsv, graph):
    colours_dict = {}
    with open(colours_tsv, "r") as f:
        #header
        f.readline()
        for line in f:
            line_list = (line.strip()).split()
            node_id = line_list[0]
            colours_dict[node_id] = line_list[1:]
    nx.set_node_attributes(graph._graph, 'colours', colours_dict)

def add_edges_to_node_attributes(graph, colours=False):
    node_pos_dict = {}
    node_neg_dict = {}
    for item in graph._graph.adjacency_iter():
        node_id, adjacency_info = item
        node_pos_dict[node_id] = {}
        node_neg_dict[node_id] = {}
        for sink_node_id, sink_node_info in adjacency_info.items():
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                #if colours are true, check that any one of the sink nodes colours matches at least one of those in the source node
                if colours == True:
                    if any(graph._graph.node[virtual_edge_info['to_node']]['colours'][i] == graph._graph.node[node_id]['colours'][i] for i in range(0, len(graph._graph.node[node_id]['colours']))):
                        if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+':
                            node_pos_dict[node_id][virtual_edge_info['to_node']] = (virtual_edge_info['to_orn'], graph._graph.node[virtual_edge_info['to_node']]['colours'])
                        elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-':
                            node_neg_dict[node_id][virtual_edge_info['to_node']] = (virtual_edge_info['to_orn'], graph._graph.node[virtual_edge_info['to_node']]['colours'])
                #if not, just add the node and it's direction
                else:
                    if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+':
                        node_pos_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                    elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-':
                        node_neg_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
    nx.set_node_attributes(graph._graph, '+', node_pos_dict)
    nx.set_node_attributes(graph._graph, '-', node_neg_dict)

#generate graph using colours or no colours file
def generate_graph(GFA, colours_file=None):
    #generate graph
    graph = pygfa.gfa.GFA.from_file(GFA)

    #add colours and node edge attributes
    if colours_file != None:
        add_colours(colours_file, graph)
        add_edges_to_node_attributes(graph, colours=True)
    else:
        add_edges_to_node_attributes(graph, colours=False)

    return graph


class Path:
    def __init__(self, GFA, nodes, ksize, codon1=None, codon2=None, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False, colours=False):
        #initialising entries for path, where nodes is a list
        self.nodes = [nodes[0]]

        #initialise orientation; absori being orientation of first node in whole path, relori being orientation of last node added
        self.absori = startdir
        self.relori = startdir

        #create edge list for beginning node
        self.edges = self.update_edges(self.nodes[0], GFA, colours)

        #add sequences for detailed nodes if they are connected to eachother. If colours is true, checks to see if next node contains at least one
        #of the same colours in the original start node.
        if self.absori == "-":
            self.seq = str(Seq(GFA._graph.node[self.nodes[0]]['sequence']).reverse_complement())
        else:
            self.seq = GFA._graph.node[self.nodes[0]]['sequence']

        #set variable for sequence length of current unitig and number of codons generated before adding new nodes
        self.len = len(self.seq)
        self.prev_codon_len = (self.len - 3 + 1)

        if len(nodes) == 1:
            pass
        else:
            for index, i in enumerate(nodes[1:]):
                edge_list = self.edges
                for edge in edge_list:
                    source, source_dir, target, target_dir = edge
                    if source == str(nodes[index]) and target == str(i):
                        self.prev_codon_len = (self.len - 3 + 1)
                        self.seq = self.merge_path(self.seq, target, target_dir, ksize, GFA)
                        self.len = len(self.seq)
                        self.nodes.append(target)
                        self.relori = target_dir
                        self.edges = self.update_edges(self.nodes[-1], GFA, colours)
                    else:
                        pass

        #generate codon k-mers for path
        #self.codons = self.update_kmers(self.seq)

        #create booleans to check if all paths are complete for iterative algorithm, update if codons are present
        self.frame1_complete = frame1_complete
        self.frame2_complete = frame2_complete
        self.frame3_complete = frame3_complete

        #check if k-mers are present depending on absolute orientiation of merged unitig based on first unitig, reversing ordering of codon indexing.
        if codon1 != None and codon2 != None:
            if frame1_complete == False:
                self.check_frame(codon1, codon2, 1)

            if frame2_complete == False:
                self.check_frame(codon1, codon2, 2)

            if frame3_complete == False:
                self.check_frame(codon1, codon2, 3)

        #check if all paths are complete
        self.all_frames_complete = False
        self.check_all_complete()

    # class method to genereate edges from input node
    def update_edges(self, node, GFA, colours):
        edge_list = []
        for sink_nodeid, sink_node_info in GFA._graph.edge[str(node)].items():
            # iterate through node info, picking out one of either virtual edges where the node and orientation of interest is the 'from_node'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(node) and virtual_edge_info['from_orn'] == self.relori:
                    # if colours is false, add the node regardless of colour
                    if colours == False:
                        node_tuple = (
                            virtual_edge_info['from_node'], virtual_edge_info['from_orn'],
                            virtual_edge_info['to_node'],
                            virtual_edge_info['to_orn'])
                        edge_list.append(node_tuple)
                    # if colours is true, check that at least one colour present in the original start node is present in the new node.
                    else:
                        colours_equal = False
                        for index, item in enumerate(GFA._graph.node[self.nodes[0]]['colours']):
                            if item == '1' and GFA._graph.node[virtual_edge_info['to_node']]['colours'][index] == '1':
                                colours_equal = True
                        if colours_equal == True:
                            node_tuple = (
                                virtual_edge_info['from_node'], virtual_edge_info['from_orn'],
                                virtual_edge_info['to_node'],
                                virtual_edge_info['to_orn'])
                            edge_list.append(node_tuple)
        return edge_list

    #check designated frame to see if it is complete. Also check if codon1 and codon2 are reverse due to reverse complementation
    def check_frame(self, codon1, codon2, frame):
        modulus = frame - 1
        indices1 = [m.start() for m in re.finditer(codon1, self.seq)]
        indices2 = [m.start() for m in re.finditer(codon2, self.seq)]
        codon1_frame = []
        codon2_frame = []
        for index in indices1:
            if index % 3 == modulus or index == modulus:
                codon1_frame.append(index)
        for index in indices2:
            if index % 3 == modulus or index == modulus:
                codon2_frame.append(index)
        #need to check whether any of codon2 are in next unitig, as nodes of length 1 will have end codon1 which is not
        #paired with an end stop codon, meaning a ORF will be missing.
        if any(codon1_frame[0] < (a - (self.prev_codon_len - 1)) for a in codon2_frame) or not codon1_frame:
            if frame == 1:
                self.frame1_complete = True
            elif frame == 2:
                self.frame2_complete = True
            elif frame == 3:
                self.frame3_complete = True
            else:
                pass

    #check if all frames are complete
    def check_all_complete(self):
        if self.frame1_complete == True and self.frame2_complete == True and self.frame3_complete == True:
            self.all_frames_complete = True

    #merge paths between nodes, and update the orientiation of the final node
    def merge_path(self, self_seq, new_node, new_seq_dir, ksize, GFA):
        merged_path = None
        new_seq = GFA._graph.node[str(new_node)]['sequence']
        if new_seq_dir == "+":
            merged_path = self_seq + new_seq[(ksize - 1)::]
        elif new_seq_dir == "-":
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = self_seq + rc_new_seq[(ksize - 1)::]
        else:
            pass
        return merged_path

    def create_ORF(self, codon1, codon2, frame):
        modulus = frame - 1
        indices1 = [m.start() for m in re.finditer(codon1, self.seq)]
        indices2 = [m.start() for m in re.finditer(codon2, self.seq)]
        codon1_frame = []
        codon2_frame = []
        for index in indices1:
            if index % 3 == modulus or index == modulus:
                codon1_frame.append(index)
        for index in indices2:
            if index % 3 == modulus or index == modulus:
                codon2_frame.append(index)

        # iterate through start and stop codon lists, starting at first start codon after first stop codon and pairing that with next following stop codon etc.
        ORF_indices = []
        if codon2_frame:
            last_codon2 = codon2_frame[0]
            for codon_index_1 in codon1_frame:
                if codon_index_1 > last_codon2:
                    for codon_index_2 in codon2_frame:
                        if codon_index_2 > codon_index_1:
                            ORF_indices.append((codon_index_1, codon_index_2))
                            last_codon2 = codon_index_2
                            break
        ORF_list = []
        for item in ORF_indices:
            start, stop = item
            ORF = self.seq[start] + self.seq[start + 1] + self.seq[start + 2]
            for i in range(start + 3, stop + 3):
                ORF = ORF + self.seq[i]
            ORF_list.append(ORF)
        return ORF_list


#recursive algorithm to generate strings and nodes with codons
def recur_paths(GFA, start_node_list, codon1, codon2, ksize, repeat, length, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False, colours=False):

    path_list = [start_node_list]
    start_path = Path(GFA, start_node_list, ksize, startdir=startdir, codon1=codon1, codon2=codon2, frame1_complete=frame1_complete, frame2_complete=frame2_complete, frame3_complete=frame3_complete, colours=colours)

    if start_path.all_frames_complete == False:
        for neighbour in start_path.edges:
            source, source_dir, target, target_dir = neighbour
            if repeat == False and target in start_path.nodes:
                pass
            else:
                path = start_node_list + [target]
                if start_path.len <= length:
                    for iteration in recur_paths(GFA, path, codon1, codon2, ksize, repeat, length, startdir=startdir, frame1_complete=start_path.frame1_complete, frame2_complete=start_path.frame2_complete, frame3_complete=start_path.frame3_complete, colours=colours):
                        path_list.append(iteration)
        del path_list[0]
    else:
        pass
    return path_list

#run recur_paths for nodes within a list
def run_recur_paths(GFA, codon1, codon2, ksize, repeat, startdir="+", length=float('inf'), colours=False):
    all_ORF_paths = {}

    #search for reverse complement of codon1 if startdir is negative, else search for codon1
    if startdir == "-":
        rc_codon1 = str(Seq(codon1).reverse_complement())
        codon1_nodes = GFA.search(lambda x: rc_codon1 in x['sequence'], limit_type=gfa.Element.NODE)
    else:
        codon1_nodes = GFA.search(lambda x: codon1 in x['sequence'], limit_type=gfa.Element.NODE)

    #run recur_paths for each item in codon1_nodes
    for node in codon1_nodes:
        print("Computing node: {}".format(node))
        node_list = [node]
        all_ORF_paths[node] = []
        all_ORF_paths[node] = recur_paths(GFA, node_list, codon1, codon2, ksize, repeat, length, startdir=startdir, frame1_complete=False, frame2_complete=False, frame3_complete=False, colours=colours)
        print("Completed node: {}".format(node))

    return all_ORF_paths

def ORF_generation(GFA, stop_codon, start_codon, ksize, repeat, length=float('inf'), colours=False):
    all_ORF_paths = {}

    #search for all nodes with stop codon with positive and negative strandedness
    #rc_codon1 = str(Seq(stop_codon).reverse_complement())
    #stop_nodes_neg = GFA.search(lambda x: rc_codon1 in x['sequence'], limit_type=gfa.Element.NODE)
    #stop_nodes_pos = GFA.search(lambda x: stop_codon in x['sequence'], limit_type=gfa.Element.NODE)

    stop_nodes_pos = ['2']
    stop_nodes_neg = ['2']


    #run recur_paths for each stop codon detected, generating ORFs from node list
    for node in stop_nodes_pos:
        print("Computing node (pos): {}".format(node))
        node_list = [node]

        #intialise all frame dictionaries
        all_ORF_paths[node] = {}
        all_ORF_paths[node]['+'] = {}
        all_ORF_paths[node]['-'] = {}

        all_ORF_paths[node]['+'][1] = []
        all_ORF_paths[node]['+'][2] = []
        all_ORF_paths[node]['+'][3] = []

        all_ORF_paths[node]['-'][1] = []
        all_ORF_paths[node]['-'][2] = []
        all_ORF_paths[node]['-'][3] = []

        #calculate positive strand paths
        stop_to_stop_paths = recur_paths(GFA, node_list, stop_codon, stop_codon, ksize, repeat, length,
                                         startdir='+', frame1_complete=False, frame2_complete=False,
                                         frame3_complete=False, colours=colours)
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, codon1=None, codon2=None, startdir="+", frame1_complete=True,
                        frame2_complete=True, frame3_complete=True, colours=False)

            #search for ORFs using create_ORF class method
            for frame in range(1,4):
                all_ORF_paths[node]['+'][frame].extend(path.create_ORF(start_codon, stop_codon, frame))

    for node in stop_nodes_neg:
        print("Computing node (neg): {}".format(node))
        node_list = [node]

        # intialise all frame dictionaries if not present in positive strand search
        if node not in all_ORF_paths:
            all_ORF_paths[node] = {}
            all_ORF_paths[node]['+'] = {}
            all_ORF_paths[node]['-'] = {}

            all_ORF_paths[node]['+'][1] = []
            all_ORF_paths[node]['+'][2] = []
            all_ORF_paths[node]['+'][3] = []

            all_ORF_paths[node]['-'][1] = []
            all_ORF_paths[node]['-'][2] = []
            all_ORF_paths[node]['-'][3] = []

        # calculate negative strand stop-stop paths
        stop_to_stop_paths = recur_paths(GFA, node_list, stop_codon, stop_codon, ksize, repeat, length,
                                         startdir='-', frame1_complete=False, frame2_complete=False,
                                         frame3_complete=False, colours=colours)
        #calculate negative strand ORFs
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, codon1=None, codon2=None, startdir="-", frame1_complete=True,
                        frame2_complete=True, frame3_complete=True, colours=False)

            for frame in range(1, 4):
                all_ORF_paths[node]['-'][frame].extend(path.create_ORF(start_codon, stop_codon, frame))

    return all_ORF_paths

###for debugging###

if __name__ == '__main__':
    from pygfa import *
    from Bio.Seq import Seq
    import re
    graph = pygfa.gfa.GFA.from_file("test3.gfa")
    add_colours("group3_SP_capsular_gene_bifrost.tsv", graph)

    #node_list = ['424', '425', '426']
    #test_path = Path(test_graph_3, node_list, 31, "ATC", "ATC")
    #test_path.create_ORF("ATG", "ATC", 1)
    node_list = ['240', '611', '447']
    test_path = Path(graph, node_list, 31, codon1="ATC", codon2="ATC", startdir="-", frame1_complete=False,
                     frame2_complete=False, frame3_complete=False, colours=True)


    #recur_paths(graph, node_list, "ATC", "ATC", 31, False, 2000, startdir="+")