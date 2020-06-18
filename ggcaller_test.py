#imports
from pygfa import *
import networkx as nx
from Bio.Seq import Seq
import copy

#generate graph
test_graph = pygfa.gfa.GFA.from_file("test.gfa")
test_graph_3 = pygfa.gfa.GFA.from_file("test3.gfa")
group3_graph = pygfa.gfa.GFA.from_file("group3_SP_capsular_gene_bifrost.gfa")

graph = group3_graph

#length graph
len_graph = len(graph._graph.node)

#choose codon to look for
#stop_codon_neg = "ATC"
#start_codon_neg = "TAC"


#search nodes for codon
#nodes_stop = graph.search(lambda x: stop_codon_neg in x['sequence'], limit_type=gfa.Element.NODE)
#nodes_start = graph.search(lambda x: stop_codon_neg in x['sequence'], limit_type=gfa.Element.NODE)

class Path:
    def __init__(self, GFA, nodes, ksize, codon1=None, codon2=None, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False):
        #initialising entries for path, where nodes is a list
        self.nodes = [nodes[0]]

        #initialise orientation; absori being orientation of first node in whole path, relori being orientation of last node added
        self.absori = startdir
        self.relori = startdir

        #create edge list for beginning node
        self.edges = self.update_edges(self.nodes[0], GFA)

        #add sequences for detailed nodes if they are connected to eachother. If repeat = False, ignore target nodes already in path and set all frame complete to true
        if self.absori == "-":
            self.seq = str(Seq(GFA._graph.node[str(self.nodes[0])]['sequence']).reverse_complement())
        else:
            self.seq = GFA._graph.node[str(self.nodes[0])]['sequence']

        if len(nodes) == 1:
            pass
        else:
            for index, i in enumerate(nodes[1:]):
                edge_list = self.edges
                for edge in edge_list:
                    source, source_dir, target, target_dir = edge
                    if source == str(nodes[index]) and target == str(i):
                        self.seq = self.merge_path(self.seq, target, target_dir, ksize, GFA)
                        self.nodes.append(target)
                        self.relori = target_dir
                        self.edges = self.update_edges(self.nodes[-1], GFA)
                    else:
                        pass
        self.len = len(self.seq)

        #generate codon k-mers for path
        self.codons = self.update_kmers(self.seq)

        #create booleans to check if all paths are complete for iterative algorithm, update if codons are present
        self.frame1_complete = frame1_complete
        self.frame2_complete = frame2_complete
        self.frame3_complete = frame3_complete

        #check if k-mers are present depending on absolute orientiation of merged unitig based on first unitig, reversing ordering of codon indexing.
        if frame1_complete == False:
            self.check_frame(codon1, codon2, 1)

        if frame2_complete == False:
            self.check_frame(codon1, codon2, 2)

        if frame3_complete == False:
            self.check_frame(codon1, codon2, 3)

        #check if all paths are complete
        self.all_frames_complete = False
        self.check_all_complete()

    #class method to genereate edges from input node
    def update_edges(self, node, GFA):
        edge_list = []
        for sink_nodeid, sink_node_info in GFA._graph.edge[str(node)].items():
            # iterate through node info, picking out one of either virtual edges where the node and orientation of interest is the 'from_node'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(node) and virtual_edge_info['from_orn'] == self.relori:
                    node_tuple = (
                        virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'],
                        virtual_edge_info['to_orn'])
                    edge_list.append(node_tuple)
        return edge_list

    # class method to enable updating of kmers in seq from last node merged
    def update_kmers(self, seq):
        kmer_list = []
        for i in range(0, (len(seq) - 2)):
            kmer = seq[i:i + 3]
            kmer_list.append(kmer)
        return kmer_list

    #check designated frame to see if it is complete. Also check if codon1 and codon2 are reverse due to reverse complementation
    def check_frame(self, codon1, codon2, frame):
        modulus = frame - 1
        indices1 = [i for i, x in enumerate(self.codons) if x == codon1]
        indices2 = [i for i, x in enumerate(self.codons) if x == codon2]
        codon1_frame = []
        codon2_frame = []
        for index in indices1:
            if index % 3 == modulus or index == modulus:
                codon1_frame.append(index)
        for index in indices2:
            if index % 3 == modulus or index == modulus:
                codon2_frame.append(index)
        if any(a > b for a in codon2_frame for b in codon1_frame) or not codon1_frame:
            if frame == 1:
                self.frame1_complete = True
            elif frame == 2:
                self.frame2_complete = True
            elif frame == 3:
                self.frame3_complete = True
            else:
                pass

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

#recursive algorithm to generate strings and nodes with codons
def recur_paths(GFA, start_node_list, codon1, codon2, ksize, repeat, length, startdir="+", frame1_complete=False, frame2_complete=False, frame3_complete=False):

    path_list = [start_node_list]
    start_path = Path(GFA, start_node_list, ksize, startdir=startdir, codon1=codon1, codon2=codon2, frame1_complete=frame1_complete, frame2_complete=frame2_complete, frame3_complete=frame3_complete)

    if start_path.all_frames_complete == False:
        for neighbour in start_path.edges:
            source, source_dir, target, target_dir = neighbour
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


###for debugging###

if __name__ == '__main__':
    from pygfa import *
    from Bio.Seq import Seq
    graph = pygfa.gfa.GFA.from_file("group3_SP_capsular_gene_bifrost.gfa")

    node_list = ['1']
    recur_paths(graph, node_list, "ATC", "ATC", 31, False, 2000, startdir="+")