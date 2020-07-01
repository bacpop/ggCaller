#imports
from pygfa import *
import networkx as nx
import re
from Bio.Seq import Seq

#takes tsv file from Bifrost query
def add_colours(colours_tsv, GFA):
    colours_dict = {}
    with open(colours_tsv, "r") as f:
        #header
        f.readline()
        for line in f:
            line_list = (line.strip()).split()
            node_id = line_list[0]
            colours_dict[node_id] = line_list[1:]
    nx.set_node_attributes(GFA._graph, 'colours', colours_dict)

def add_edges_to_node_attributes(GFA, colours=False):
    node_pos_dict = {}
    node_neg_dict = {}
    for item in GFA._graph.adjacency_iter():
        node_id, adjacency_info = item
        node_pos_dict[node_id] = {}
        node_neg_dict[node_id] = {}
        for sink_node_id, sink_node_info in adjacency_info.items():
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                #if colours are true, check that any one of the sink nodes colours matches at least one of those in the source node and that from and to orn match
                if colours == True:
                    if any(GFA._graph.node[virtual_edge_info['to_node']]['colours'][i] == '1' and GFA._graph.node[node_id]['colours'][i] == '1' for i in range(0, len(GFA._graph.node[node_id]['colours']))):
                        if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+' and str(virtual_edge_info['to_orn']) == '+':
                            node_pos_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                        elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-' and str(virtual_edge_info['to_orn']) == '-':
                            node_neg_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                #if not, just add the node and it's direction
                else:
                    if str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '+' and str(virtual_edge_info['to_orn']) == '+':
                        node_pos_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
                    elif str(virtual_edge_info['from_node']) == str(node_id) and str(virtual_edge_info['from_orn']) == '-' and str(virtual_edge_info['to_orn']) == '-':
                        node_neg_dict[node_id][virtual_edge_info['to_node']] = virtual_edge_info['to_orn']
    nx.set_node_attributes(GFA._graph, '+', node_pos_dict)
    nx.set_node_attributes(GFA._graph, '-', node_neg_dict)

#create binary arrays for each node for codon enumeration
def enumerate_codon(GFA, codon_list, ksize):
    #initialise full and partial binary dictionaries
    full_dict = {}
    part_dict = {}

    #generate reverse complement codon list
    rev_codon_list = [str(Seq(codon).reverse_complement()) for codon in codon_list]

    #compile regular express for codon index searching
    pos_regex = re.compile((r'(' + '|'.join(codon_list) + r')'))
    neg_regex = re.compile((r'(' + '|'.join(rev_codon_list) + r')'))

    #add to each dict for each node
    for node, item in GFA._graph.node.items():
        #initialise node sequence
        seq = item['sequence']

        #enumerate codons in full sequence. For negative, indices must start from end of sequence.
        pos_indices = [m.start() for m in pos_regex.finditer(seq)]
        neg_indices = [(len(seq) - m.start() - 3) for m in neg_regex.finditer(seq)]
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
        full_array = {'+': [pos_frame1, pos_frame2, pos_frame3], '-': [neg_frame1, neg_frame2, neg_frame3]}
        full_dict[node] = full_array

        #enumerate codons in partial sequences
        part_pos_seq = seq[ksize-1:]
        part_neg_seq = seq[:-(ksize-1)]

        pos_indices = [m.start() for m in pos_regex.finditer(part_pos_seq)]
        neg_indices = [(len(part_neg_seq) - m.start() - 3) for m in neg_regex.finditer(part_neg_seq)]
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

        #create part array for three different reading frames, depending on length of previous unitig.
        part_array = {0: {'+': [pos_frame1, pos_frame2, pos_frame3], '-': [neg_frame1, neg_frame2, neg_frame3]},
                      1: {'+': [pos_frame3, pos_frame1, pos_frame2], '-': [neg_frame3, neg_frame1, neg_frame2]},
                      2: {'+': [pos_frame2, pos_frame3, pos_frame1], '-': [neg_frame2, neg_frame3, neg_frame1]}}
        part_dict[node] = part_array

    #add full and part binary dicts to GFA
    nx.set_node_attributes(GFA._graph, 'full_bin', full_dict)
    nx.set_node_attributes(GFA._graph, 'part_bin', part_dict)


#generate graph using colours or no colours file
def generate_graph(GFA, ksize, codon_list=None, colours_file=None):
    #generate graph
    graph = pygfa.gfa.GFA.from_file(GFA)

    #add colours and node edge attributes
    if colours_file is not None:
        add_colours(colours_file, graph)
        add_edges_to_node_attributes(graph, colours=True)
    else:
        add_edges_to_node_attributes(graph, colours=False)

    if codon_list is not None:
        enumerate_codon(graph, codon_list, ksize)

    return graph


class Path:
    def __init__(self, GFA, nodes, ksize, startdir="+", create_ORF=False):
        #initialising entries for path, where nodes is a list
        self.nodes = [nodes[0]]

        #initialise orientation; absori being orientation of first node in whole path, relori being orientation of last node added
        self.absori = startdir
        #self.relori = startdir

        #create edge list for beginning node
        self.edges = GFA._graph.node[self.nodes[0]][self.absori]

        #create source node colours object
        self.source_colour = GFA._graph.node[self.nodes[0]]['colours']

        #get node length and unitig frame modulus for appending node frame
        self.len = GFA._graph.node[self.nodes[0]]['slen']
        self.modulus = self.len % 3

        #get full binary matrix for start node
        full_binary_matrix = GFA._graph.node[self.nodes[0]]['full_bin'][self.absori][:]

        #if create_ORF is false, run standard methods to check if frames are complete
        if create_ORF == False:
            if len(nodes) == 1:
                pass
            else:
                for node in nodes[1:]:
                    try:
                        #check if nodes contains any of the same colours as sink node in path - gives recursion limit error
                        #if any(i == '1' and j == '1' for i, j in zip(self.source_colour, GFA._graph.node[node]['colours'])):

                        #update relative orientiation based on newly added node - removed as results in non-true sequences
                        #self.relori = self.edges[node]
                        #get part binary matrix of node to be added on, depending on length of previous nodes
                        part_binary_matrix = GFA._graph.node[node]['part_bin'][self.modulus][self.absori]

                        #conduct binary matrix subtraction to determine whether frames are complete
                        for i in range(0, 3):
                            if full_binary_matrix[i] == 1 and part_binary_matrix[i] == 1:
                                full_binary_matrix[i] = 0

                        #update path length, edges list, modulus and nodes list
                        self.len += (GFA._graph.node[node]['slen'] - (ksize - 1))
                        self.modulus = self.len % 3
                        self.edges = GFA._graph.node[node][self.absori]
                        self.nodes.append(node)
                    except KeyError:
                        pass


            #check if all paths are complete
            self.all_frames_complete = False
            if sum(full_binary_matrix) == 0:
                self.all_frames_complete = True

        #run original merge path method to generate ORF if create_ORF is true, ignoring if frames are complete
        else:
            self.path_colour = GFA._graph.node[self.nodes[0]]['colours'][:]
            if self.absori == "-":
                self.seq = str(Seq(GFA._graph.node[self.nodes[0]]['sequence']).reverse_complement())
            else:
                self.seq = GFA._graph.node[self.nodes[0]]['sequence']
            if len(nodes) == 1:
                pass
            else:
                for node in nodes[1:]:
                    try:
                        #target_dir = self.edges[node]
                        self.seq = self.merge_path(self.seq, node, self.absori, ksize, GFA)
                        self.nodes.append(node)
                        #self.relori = target_dir
                        self.edges = GFA._graph.node[self.nodes[-1]][self.absori]
                        edge_colour = GFA._graph.node[node]['colours']
                        #calculate shared colours through entire path length
                        for i in range(0, len(self.path_colour)):
                            if self.path_colour[i] == '1' and edge_colour[i] == '0':
                                self.path_colour[i] = '0'
                    except KeyError:
                        pass

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

    def create_ORF(self, start_codon_list, stop_codon_list, frame):
        #initialise frames and start/stop regular expressions
        modulus = frame - 1
        start_regex = re.compile((r'(' + '|'.join(start_codon_list) + r')'))
        stop_regex = re.compile((r'(' + '|'.join(stop_codon_list) + r')'))

        #generate indices for all start and stop codons
        start_indices = [m.start() for m in start_regex.finditer(self.seq)]
        stop_indices = [m.start() for m in stop_regex.finditer(self.seq)]

        #extract only indices that are in the same frame
        codon1_frame = [index for index in start_indices if index % 3 == modulus or index == modulus]
        codon2_frame = [index for index in stop_indices if index % 3 == modulus or index == modulus]

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
def recur_paths(GFA, start_node_list, ksize, repeat, length, startdir="+"):

    path_list = [start_node_list]
    start_path = Path(GFA, start_node_list, ksize, startdir=startdir, create_ORF=False)

    if start_path.all_frames_complete == False:
        for target in start_path.edges:
            if repeat == False and target in start_path.nodes:
                pass
            else:
                path = start_node_list + [target]
                if start_path.len <= length:
                    for iteration in recur_paths(GFA, path, ksize, repeat, length, startdir=startdir):
                        path_list.append(iteration)
        del path_list[0]
    else:
        pass
    return path_list

def ORF_generation(GFA, stop_codon_list, start_codon_list, ksize, repeat, path_freq=0, length=float('inf')):
    all_ORF_paths = {}
    all_ORF_paths['+'] = {}
    all_ORF_paths['-'] = {}

    #generate list of stop codon list
    rev_stop_codon_list = [str(Seq(codon).reverse_complement()) for codon in stop_codon_list]
    # search for all nodes with stop codon with positive and negative strandedness, add to set
    stop_nodes_pos = set()
    stop_nodes_neg = set()

    for codon in stop_codon_list:
        stop_nodes_pos.update(tuple(GFA.search(lambda x: codon in x['sequence'], limit_type=gfa.Element.NODE)))
    for codon in rev_stop_codon_list:
        stop_nodes_neg.update(tuple(GFA.search(lambda x: codon in x['sequence'], limit_type=gfa.Element.NODE)))

    #calculate length of lists for standard-out
    len_pos_list = len(stop_nodes_pos)
    len_neg_list = len(stop_nodes_neg)

    #run recur_paths for each stop codon detected in positive list, generating ORFs from node list
    count = 0
    for node in stop_nodes_pos:
        count += 1
        print("Computing positive node: {} / {}".format(count, len_pos_list))
        node_list = [node]

        #calculate positive strand paths
        stop_to_stop_paths = recur_paths(GFA, node_list, ksize, repeat, length, startdir="+")
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, startdir="+", create_ORF=True)

            #check to see if path contains number of colours greater than frequency specified
            path_colour_freq = (path.path_colour.count('1'))/len(path.path_colour)
            if path_colour_freq >= path_freq and path_colour_freq > 0:
                #search for ORFs using create_ORF class method
                for frame in range(1, 4):
                    #search for each start codon and stop codon combination
                    for start_codon in start_codon_list:
                        start_codon = [start_codon]
                        ORF_list = path.create_ORF(start_codon, stop_codon_list, frame)
                        for ORF in ORF_list:
                            if ORF not in all_ORF_paths['+']:
                                all_ORF_paths['+'][ORF] = [path.path_colour]
                            else:
                                all_ORF_paths['+'][ORF].append(path.path_colour)

    # run recur_paths for each stop codon detected in negative list, generating ORFs from node list
    count = 0
    for node in stop_nodes_neg:
        count += 1
        print("Computing negative node: {} / {}".format(count, len_neg_list))
        node_list = [node]

        # calculate negative strand paths
        stop_to_stop_paths = recur_paths(GFA, node_list, ksize, repeat, length, startdir="-")
        for node_path in stop_to_stop_paths:
            path = Path(GFA, node_path, ksize, startdir="-", create_ORF=True)

            # check to see if path contains number of colours greater than frequency specified
            path_colour_freq = (path.path_colour.count('1')) / len(path.path_colour)
            if path_colour_freq >= path_freq and path_colour_freq > 0:
                # search for ORFs using create_ORF class method
                for frame in range(1, 4):
                    # search for each start codon and stop codon combination
                    for start_codon in start_codon_list:
                        start_codon = [start_codon]
                        ORF_list = path.create_ORF(start_codon, stop_codon_list, frame)
                        for ORF in ORF_list:
                            if ORF not in all_ORF_paths['-']:
                                all_ORF_paths['-'][ORF] = [path.path_colour]
                            else:
                                all_ORF_paths['-'][ORF].append(path.path_colour)

    return all_ORF_paths

###for debugging###

if __name__ == '__main__':
    import sys
    from pygfa import *
    from Bio.Seq import Seq
    import re
    import networkx

    #for debugging
    #stop_codon_list = ["TAA", "TGA", "TAG"]
    #start_codon_list = ["ATG", "GTG", "TTG"]
    #graph = generate_graph("group3_SP_capsular_gene_bifrost.gfa", 31, stop_codon_list, "group3_SP_capsular_gene_bifrost.tsv")
    #ORF_generation(graph, stop_codon_list, start_codon_list, 31, False, length=2000)


    graph_file = sys.argv[1]
    tsv_file = sys.argv[2]
    ksize = int(sys.argv[3])
    path_length = int(sys.argv[4])
    ORF_length = int(sys.argv[5])
    output = sys.argv[6]

    stop_codon_list = ["TAA", "TGA", "TAG"]
    start_codon_list = ["ATG", "GTG", "TTG"]
    graph = generate_graph(graph_file, ksize, stop_codon_list, tsv_file)

    ORF_output = ORF_generation(graph, stop_codon_list, start_codon_list, ksize, False, length=path_length)

    with open(output, "w") as f:
        count = 1
        for key, item in ORF_output['+'].items():
            if len(key) >= ORF_length:
                f.write(">Gene_ID: " + str(count) + " Strand: +" + "\n" + str(key) + "\n")
                count += 1
        for key, item in ORF_output['-'].items():
            if len(key) >= ORF_length:
                f.write(">Gene_ID: " + str(count) + " Strand: -" + "\n" + str(key) + "\n")
                count += 1