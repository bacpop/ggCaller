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

def add_edges_to_node_attributes(GFA, colours=False):
    node_pos_dict = {}
    node_neg_dict = {}
    for item in GFA._graph.adjacency_iter():
        node_id, adjacency_info = item
        node_pos_dict[node_id] = {}
        node_neg_dict[node_id] = {}
        for sink_node_id, sink_node_info in adjacency_info.items():
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                #if colours are true, check that any one of the sink nodes colours matches at least one of those in the source node
                if colours == True:
                    if any(GFA._graph.node[virtual_edge_info['to_node']]['colours'][i] == GFA._graph.node[node_id]['colours'][i] for i in range(0, len(GFA._graph.node[node_id]['colours']))):
                        if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+':
                            node_pos_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                        elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-':
                            node_neg_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                #if not, just add the node and it's direction
                else:
                    if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+':
                        node_pos_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                    elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-':
                        node_neg_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
    nx.set_node_attributes(GFA._graph, '+', node_pos_dict)
    nx.set_node_attributes(GFA._graph, '-', node_neg_dict)

#create binary arrays for each node for codon enumeration
def enumerate_codon(GFA, codon, ksize):
    #initialise full and partial binary dictionaries
    full_dict = {}
    part_dict = {}

    #generate reverse complement codon
    rev_codon = str(Seq(codon).reverse_complement())

    #add to each dict for each node
    for node, item in GFA._graph.node.items():
        #initialise node sequence
        seq = item['sequence']

        #enumerate codons in full sequence. For negative, indices must start from end of sequence.
        pos_indices = [m.start() for m in re.finditer(codon, seq)]
        neg_indices = [(len(seq) - m.start() - 3) for m in re.finditer(rev_codon, seq)]
        pos_frame1 = 0
        pos_frame2 = 0
        pos_frame3 = 0
        neg_frame1 = 0
        neg_frame2 = 0
        neg_frame3 = 0
        for index in pos_indices:
            if index % 3 == 0 or index == 0:
                pos_frame1 = 1
            elif index % 3 == 1 or index == 1:
                pos_frame2 = 1
            elif index % 3 == 2 or index == 2:
                pos_frame3 = 1
            else:
                pass
        for index in neg_indices:
            if index % 3 == 0 or index == 0:
                neg_frame1 = 1
            elif index % 3 == 1 or index == 1:
                neg_frame2 = 1
            elif index % 3 == 2 or index == 2:
                neg_frame3 = 1
            else:
                pass

        #create full array, add it to dictionary
        full_array = [pos_frame1, pos_frame2, pos_frame3, neg_frame1, neg_frame2, neg_frame3]
        full_dict[node] = full_array

        #enumerate codons in partial sequences
        part_pos_seq = seq[ksize-1:]
        part_neg_seq = seq[:-(ksize-1)]

        pos_indices = [m.start() for m in re.finditer(codon, part_pos_seq)]
        neg_indices = [(len(part_neg_seq) - m.start() - 3) for m in re.finditer(rev_codon, part_neg_seq)]
        pos_frame1 = 0
        pos_frame2 = 0
        pos_frame3 = 0
        neg_frame1 = 0
        neg_frame2 = 0
        neg_frame3 = 0
        for index in pos_indices:
            if index % 3 == 0 or index == 0:
                pos_frame1 = 1
            elif index % 3 == 1 or index == 1:
                pos_frame2 = 1
            elif index % 3 == 2 or index == 2:
                pos_frame3 = 1
            else:
                pass
        for index in neg_indices:
            if index % 3 == 0 or index == 0:
                neg_frame1 = 1
            elif index % 3 == 1 or index == 1:
                neg_frame2 = 1
            elif index % 3 == 2 or index == 2:
                neg_frame3 = 1
            else:
                pass

        part_array = [pos_frame1, pos_frame2, pos_frame3, neg_frame1, neg_frame2, neg_frame3]
        part_dict[node] = part_array

    #add full and part binary dicts to GFA
    nx.set_node_attributes(GFA._graph, 'full_bin', full_dict)
    nx.set_node_attributes(GFA._graph, 'part_bin', part_dict)


#generate graph using colours or no colours file
def generate_graph(GFA, ksize, codon=None, colours_file=None):
    #generate graph
    graph = pygfa.gfa.GFA.from_file(GFA)

    #add colours and node edge attributes
    if colours_file is not None:
        add_colours(colours_file, graph)
        add_edges_to_node_attributes(graph, colours=True)
    else:
        add_edges_to_node_attributes(graph, colours=False)

    if codon is not None:
        enumerate_codon(graph, codon, ksize)

    return graph


class Path:
    def __init__(self, GFA, nodes, ksize, codon1=None, codon2=None, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False, create_ORF=False):
        #initialising entries for path, where nodes is a list
        self.nodes = [nodes[0]]

        #initialise orientation; absori being orientation of first node in whole path, relori being orientation of last node added
        self.absori = startdir
        self.relori = startdir

        #create edge list for beginning node
        self.edges = GFA._graph.node[self.nodes[0]][self.absori]

        #if create_ORF is false, run standard methods to check if frames are complete
        if create_ORF == False:
            #add sequences for detailed nodes if they are connected to eachother. If colours is true, checks to see if next node contains at least one
            #of the same colours in the original start node.
            if self.absori == "-":
                self.iter_seq = str(Seq(GFA._graph.node[self.nodes[0]]['sequence']).reverse_complement())
            else:
                self.iter_seq = GFA._graph.node[self.nodes[0]]['sequence']

            #set variable for sequence length of current unitig and number of codons generated before adding new nodes
            self.len = len(self.iter_seq)
            self.prev_codon_len = (self.len - 3 + 1)

            #create indices for codons in each frame of self.prev_seq
            indices1 = [m.start() for m in re.finditer(codon1, self.iter_seq)]
            indices2 = [m.start() for m in re.finditer(codon2, self.iter_seq)]

            self.codon1_frame1 = [index for index in indices1 if index % 3 == 0 or index == 0]
            self.codon2_frame1 = [index for index in indices2 if index % 3 == 0 or index == 0]
            self.codon1_frame2 = [index for index in indices1 if index % 3 == 1 or index == 1]
            self.codon2_frame2 = [index for index in indices2 if index % 3 == 1 or index == 1]
            self.codon1_frame3 = [index for index in indices1 if index % 3 == 2 or index == 2]
            self.codon2_frame3 = [index for index in indices2 if index % 3 == 2 or index == 2]

            #iterate through nodes to be added, adding them in their correct orientations and updating self.edges
            if len(nodes) == 1:
                pass
            else:
                for node in nodes[1:]:
                    try:
                        target_dir = self.edges[node]
                        self.prev_codon_len = (self.len - 3 + 1)
                        #generate indices for codon1 and codon 2
                        self.iterate_sequence(node, target_dir, ksize, codon1, codon2, GFA)
                        self.nodes.append(node)
                        self.relori = target_dir
                        self.edges = GFA._graph.node[self.nodes[-1]][self.relori]
                    except KeyError:
                        pass


            #create booleans to check if all paths are complete for iterative algorithm, update if codons are present
            self.frame1_complete = frame1_complete
            self.frame2_complete = frame2_complete
            self.frame3_complete = frame3_complete

            #check if k-mers are present depending on absolute orientiation of merged unitig based on first unitig, reversing ordering of codon indexing.
            if codon1 != None and codon2 != None:
                if self.frame1_complete == False:
                    self.check_frame(1)

                if self.frame2_complete == False:
                    self.check_frame(2)

                if self.frame3_complete == False:
                    self.check_frame(3)

            #check if all paths are complete
            self.all_frames_complete = False
            self.check_all_complete()

        #run original merge path method to generate ORF if create_ORF is true, ignoring if frames are complete
        else:
            if self.absori == "-":
                self.seq = str(Seq(GFA._graph.node[self.nodes[0]]['sequence']).reverse_complement())
            else:
                self.seq = GFA._graph.node[self.nodes[0]]['sequence']
            if len(nodes) == 1:
                pass
            else:
                for node in nodes[1:]:
                    try:
                        target_dir = self.edges[node]
                        self.seq = self.merge_path(self.seq, node, target_dir, ksize, GFA)
                        self.nodes.append(node)
                        self.relori = target_dir
                        self.edges = GFA._graph.node[self.nodes[-1]][self.relori]
                    except KeyError:
                        pass

    def iterate_sequence(self, node, target_dir, ksize, codon1, codon2, GFA):
        if target_dir == "-":
            self.iter_seq = str(Seq(GFA._graph.node[node]['sequence']).reverse_complement())
        else:
            self.iter_seq = GFA._graph.node[node]['sequence']

        # generate indices of all codons in iter_seq
        indices1 = [m.start() for m in re.finditer(codon1, self.iter_seq)]
        indices2 = [m.start() for m in re.finditer(codon2, self.iter_seq)]

        #create a new list for each index in new sequence, ensuring that its index is added to the self.len and that it is not repeated in the previous list
        codon1_frame1 = [(index + self.len - (ksize - 1)) for index in indices1 if ((index + self.len - (ksize - 1)) % 3 == 0 or index == 0) and index > (ksize - 4)]
        codon2_frame1 = [(index + self.len - (ksize - 1)) for index in indices2 if ((index + self.len - (ksize - 1)) % 3 == 0 or index == 0) and index > (ksize - 4)]
        codon1_frame2 = [(index + self.len - (ksize - 1)) for index in indices1 if ((index + self.len - (ksize - 1)) % 3 == 1 or index == 1) and index > (ksize - 4)]
        codon2_frame2 = [(index + self.len - (ksize - 1)) for index in indices2 if ((index + self.len - (ksize - 1)) % 3 == 1 or index == 1) and index > (ksize - 4)]
        codon1_frame3 = [(index + self.len - (ksize - 1)) for index in indices1 if ((index + self.len - (ksize - 1)) % 3 == 2 or index == 2) and index > (ksize - 4)]
        codon2_frame3 = [(index + self.len - (ksize - 1)) for index in indices2 if ((index + self.len - (ksize - 1)) % 3 == 2 or index == 2) and index > (ksize - 4)]

        #append new indices to each codon1/2 frame variables
        self.codon1_frame1 += codon1_frame1
        self.codon2_frame1 += codon2_frame1
        self.codon1_frame2 += codon1_frame2
        self.codon2_frame2 += codon2_frame2
        self.codon1_frame3 += codon1_frame3
        self.codon2_frame3 += codon2_frame3

        #update self.len variable to include length of new unitig, minus overlap
        self.len = self.len + (len(self.iter_seq) - (ksize - 1))

        #check designated frame to see if it is complete. Also check if codon1 and codon2 are reverse due to reverse complementation
    def check_frame(self, frame):
        if frame == 1:
            codon1_frame = self.codon1_frame1
            codon2_frame = self.codon2_frame1
        elif frame == 2:
            codon1_frame = self.codon1_frame2
            codon2_frame = self.codon2_frame2
        elif frame == 3:
            codon1_frame = self.codon1_frame3
            codon2_frame = self.codon2_frame3
        else:
            pass

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

    # class method to genereate edges from input node
    def update_edges_old(self, node, GFA, colours):
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

        codon1_frame = [index for index in indices1 if index % 3 == modulus or index == modulus]
        codon2_frame = [index for index in indices2 if index % 3 == modulus or index == modulus]

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
            ORF = self.seq[start:(stop + 3)]
            ORF_list.append(ORF)
        return ORF_list


#recursive algorithm to generate strings and nodes with codons
def recur_paths(GFA, start_node_list, codon1, codon2, ksize, repeat, length, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False):

    path_list = [start_node_list]
    start_path = Path(GFA, start_node_list, ksize, startdir=startdir, codon1=codon1, codon2=codon2, frame1_complete=frame1_complete, frame2_complete=frame2_complete, frame3_complete=frame3_complete, create_ORF=False)

    if start_path.all_frames_complete == False:
        for target in start_path.edges:
            if repeat == False and target in start_path.nodes:
                pass
            else:
                path = start_node_list + [target]
                if start_path.len <= length:
                    for iteration in recur_paths(GFA, path, codon1, codon2, ksize, repeat, length, startdir=startdir, frame1_complete=start_path.frame1_complete, frame2_complete=start_path.frame2_complete, frame3_complete=start_path.frame3_complete):
                        path_list.append(iteration)
        del path_list[0]
    else:
        pass
    return path_list

#run recur_paths for nodes within a list
def run_recur_paths(GFA, codon1, codon2, ksize, repeat, startdir="+", length=float('inf')):
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
        all_ORF_paths[node] = recur_paths(GFA, node_list, codon1, codon2, ksize, repeat, length, startdir=startdir, frame1_complete=False, frame2_complete=False, frame3_complete=False)
        print("Completed node: {}".format(node))

    return all_ORF_paths

def ORF_generation(GFA, stop_codon, start_codon, ksize, repeat, length=float('inf')):
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
                                         frame3_complete=False)
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, codon1=None, codon2=None, startdir="+", frame1_complete=True,
                        frame2_complete=True, frame3_complete=True, create_ORF=True)

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
                                         frame3_complete=False)
        #calculate negative strand ORFs
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, codon1=None, codon2=None, startdir="-", frame1_complete=True,
                        frame2_complete=True, frame3_complete=True, create_ORF=True)

            for frame in range(1, 4):
                all_ORF_paths[node]['-'][frame].extend(path.create_ORF(start_codon, stop_codon, frame))

    return all_ORF_paths

###for debugging###

if __name__ == '__main__':
    from pygfa import *
    from Bio.Seq import Seq
    import re

    graph = generate_graph("test3.gfa", 31, "ATC", "group3_SP_capsular_gene_bifrost.tsv")
