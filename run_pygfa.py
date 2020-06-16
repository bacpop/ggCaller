#imports
from pygfa import *
import networkx as nx
from Bio.Seq import Seq
import copy

#generate graph
graph = pygfa.gfa.GFA.from_file("test.gfa")

#length graph
len_graph = len(graph._graph.node)

#choose codon to look for
stop_codon_neg = "ACT"
start_codon_neg = "TAC"


#search nodes for codon
nodes_stop = graph.search(lambda x: stop_codon_neg in x['sequence'], limit_type=gfa.Element.NODE)
nodes_start = graph.search(lambda x: stop_codon_neg in x['sequence'], limit_type=gfa.Element.NODE)

#method to get other nodes connected to a node, and the orientiation of the edges (note this is bidirected graph, so edges can be traversed in either way via reverse complement)
def get_edges_from_node(GFA, node):
    edge_list = []
    for sink_nodeid, sink_node_info in GFA._graph.edge[str(node)].items():
        #iterate through node info, picking out one of either virtual edges where the node of interest is the 'from_orn'
        for virtual_edgeid, virtual_edge_info in sink_node_info.items():
            if virtual_edge_info['from_node'] == str(node):
                node_tuple = (virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'], virtual_edge_info['to_orn'])
                edge_list.append(node_tuple)
    return edge_list

#use nx dijkstra_path to compute paths between nodes
def run_dijkstra(in_graph, node_list1, node_list2):
    dij_paths = {}
    for i in node_list1:
        visited = set()
        dij_paths[i] = []
        for j in node_list2:
            if i == j:
                pass
            else:
                if j not in visited:
                    list_paths = nx.dijkstra_path(in_graph, i, j)
                    dij_paths[i].append(list_paths)
                    visited.update(list_paths)
    return dij_paths

def merge_unitigs(unitig1, dir1, unitig2, dir2, ksize):
    merged_unitig = None
    if dir1 == "+" and dir2 == "+":
        merged_unitig = unitig1 + unitig2[(ksize - 1)::]
    elif dir1 == "+" and dir2 == "-":
        seq2 = Seq(unitig2)
        rc_seq2 = str(seq2.reverse_complement())
        merged_unitig = unitig1 + rc_seq2[(ksize - 1)::]
    elif dir1 == "-" and dir2 == "+":
        seq1 = Seq(unitig1)
        rc_seq1 = str(seq1.reverse_complement())
        merged_unitig = unitig2 + rc_seq1[(ksize - 1)::]
    elif dir1 == "-" and dir2 == "-":
        seq1 = Seq(unitig1)
        rc_seq1 = str(seq1.reverse_complement())
        seq2 = Seq(unitig2)
        rc_seq2 = str(seq2.reverse_complement())
        merged_unitig = rc_seq2 + rc_seq1[(ksize - 1)::]
    else:
        pass
    return merged_unitig

def find_codon_indices(codon, unitig):
    import re
    codon_indices = []
    for seq in re.finditer(codon, unitig):
        codon_indices.append(seq.start())
    return codon_indices

def generate_kmers(unitig, k):
    kmer_list = []
    for i in range(0, (len(unitig) - (k-1))):
        kmer = unitig[i:i+k]
        kmer_list.append(kmer)
    return kmer_list

#create a class called path, which enables addition of connected unitigs
class Path:
    def __init__(self, id, node, graph):
        #initialising entries for path
        self.id = id
        self.nodes = [str(node)]
        self.seq = graph._graph.node[str(node)]['sequence']
        self.len = len(graph._graph.node[str(node)]['sequence'])
        self.ori = "+"

        #generate codon k-mers for path
        kmer_list = []
        for i in range(0, (len(self.seq) - 3)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

        #create booleans to check if all paths are complete for iterative algorithm
        self.frame1_complete = False
        self.frame2_complete = False
        self.frame3_complete = False
        self.all_frames_complete = False

        #create edge list for end node
        edge_list = []
        for sink_nodeid, sink_node_info in graph._graph.edge[str(self.nodes[-1])].items():
            # iterate through node info, picking out one of either virtual edges where the node of interest is the 'from_orn'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(self.nodes[-1]):
                    node_tuple = (
                        virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'],
                        virtual_edge_info['to_orn'])
                    edge_list.append(node_tuple)
        self.edges = edge_list

    #class method to enable updating of edges from last node merged
    def update_edges(self):
        edge_list = []
        for sink_nodeid, sink_node_info in graph._graph.edge[str(self.nodes[-1])].items():
            # iterate through node info, picking out one of either virtual edges where the node of interest is the 'from_orn'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(self.nodes[-1]):
                    node_tuple = (
                        virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'],
                        virtual_edge_info['to_orn'])
                    edge_list.append(node_tuple)
        self.edges = edge_list

    # class method to enable updating of kmers in seq from last node merged
    def updates_kmers(self):
        kmer_list = []
        for i in range(0, (len(self.seq) - 3)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

    #allow ability to edit path id
    def change_id(self, name):
        self.id = name

    #allow ability to reset frame complete status
    def reset_frame(self):
        self.frame1_complete = False
        self.frame2_complete = False
        self.frame3_complete = False
        self.all_frames_complete = False

    def index_kmers(self, codon1, codon2):
        indices1 = [i for i, x in enumerate(self.codons) if x == codon1]
        indices2 = [i for i, x in enumerate(self.codons) if x == codon2]
        codon1_frame1 = []
        codon1_frame2 = []
        codon1_frame3 = []
        codon2_frame1 = []
        codon2_frame2 = []
        codon2_frame3 = []
        for index in indices1:
            if index % 3 == 0:
                codon1_frame1.append(index)
            elif index % 3 == 1:
                codon1_frame2.append(index)
            elif index % 3 == 2:
                codon1_frame3.append(index)
            else:
                pass
        for index in indices2:
            if index % 3 == 0:
                codon2_frame1.append(index)
            elif index % 3 == 1:
                codon2_frame2.append(index)
            elif index % 3 == 2:
                codon2_frame3.append(index)
            else:
                pass

        # if any of the codon2 indexes are higher than those in the codon1 for each frame, change the frames to complete
        if any(a > b for a in codon2_frame1 for b in codon1_frame1):
            self.frame1_complete = True
        if any(a > b for a in codon2_frame2 for b in codon1_frame2):
            self.frame2_complete = True
        if any(a > b for a in codon2_frame3 for b in codon1_frame3):
            self.frame3_complete = True
        #if all frames are complete, set all_frames_complete to true
        if self.frame1_complete == True and self.frame2_complete == True and self.frame3_complete == True:
            self.all_frames_complete = True

    #merge paths between nodes, and update the orientiation of the final node
    def merge_path(self, self_dir, new_node, new_seq_dir, ksize, graph):
        merged_path = None
        new_seq = graph._graph.node[str(new_node)]['sequence']
        if self_dir == "+" and new_seq_dir == "+":
            merged_path = self.seq + new_seq[(ksize - 1)::]
            self.ori = "+"
        elif self_dir == "+" and new_seq_dir == "-":
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = self.seq + rc_new_seq[(ksize - 1)::]
            self.ori = "-"
        elif self_dir == "-" and new_seq_dir == "+":
            seq1 = Seq(self.seq)
            rc_seq1 = str(seq1.reverse_complement())
            merged_path = new_seq + rc_seq1[(ksize - 1)::]
            self.ori = "+"
        elif self_dir == "-" and new_seq_dir == "-":
            seq1 = Seq(self.seq)
            rc_seq1 = str(seq1.reverse_complement())
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = rc_new_seq + rc_seq1[(ksize - 1)::]
            self.ori = "-"
        else:
            pass
        self.seq = merged_path
        self.len = len(merged_path)
        self.nodes.append(new_node)
        self.update_edges()
        self.updates_kmers()

class Path2:
    def __init__(self, GFA, id, nodes, ksize, codon1 = None, codon2 = None):
        #initialising entries for path, where nodes is a list
        self.id = id
        self.nodes = [nodes[0]]
        self.ori = "+"

        #create edge list for beginning node
        self.edges = self.update_edges(self.nodes[0], GFA)

        #add sequences for detailed nodes if they are connected to eachother
        self.seq = GFA._graph.node[str(self.nodes[0])]['sequence']
        if len(nodes) == 1:
            pass
        else:
            for index, i in enumerate(nodes[1:]):
                edge_list = self.edges
                for edge in edge_list:
                    source, source_dir, target, target_dir = edge
                    if source == str(nodes[index]) and target == str(i):
                        self.seq = self.merge_path(self.seq, source_dir, target, target_dir, ksize, GFA)
                        self.nodes.append(target)
                        self.edges = self.update_edges(self.nodes[-1], GFA)
                    else:
                        pass
        self.len = len(self.seq)

        #generate codon k-mers for path
        kmer_list = []
        for i in range(0, (len(self.seq) - 3)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

        #create booleans to check if all paths are complete for iterative algorithm, update if codons are present
        self.frame1_complete = False
        self.frame2_complete = False
        self.frame3_complete = False
        self.all_frames_complete = False
        self.index_kmers(codon1, codon2)

    #class method to genereate edges from input node
    def update_edges(self, node, GFA):
        edge_list = []
        for sink_nodeid, sink_node_info in GFA._graph.edge[str(node)].items():
            # iterate through node info, picking out one of either virtual edges where the node of interest is the 'from_orn'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(node):
                    node_tuple = (
                        virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'],
                        virtual_edge_info['to_orn'])
                    edge_list.append(node_tuple)
        return edge_list

    # class method to enable updating of kmers in seq from last node merged
    def updates_kmers(self):
        kmer_list = []
        for i in range(0, (len(self.seq) - 3)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

    #allow ability to edit path id
    def change_id(self, name):
        self.id = name

    #allow ability to reset frame complete status
    def reset_frame(self):
        self.frame1_complete = False
        self.frame2_complete = False
        self.frame3_complete = False
        self.all_frames_complete = False

    def index_kmers(self, codon1, codon2):
        indices1 = [i for i, x in enumerate(self.codons) if x == codon1]
        indices2 = [i for i, x in enumerate(self.codons) if x == codon2]
        codon1_frame1 = []
        codon1_frame2 = []
        codon1_frame3 = []
        codon2_frame1 = []
        codon2_frame2 = []
        codon2_frame3 = []
        for index in indices1:
            if index % 3 == 0:
                codon1_frame1.append(index)
            elif index % 3 == 1:
                codon1_frame2.append(index)
            elif index % 3 == 2:
                codon1_frame3.append(index)
            else:
                pass
        for index in indices2:
            if index % 3 == 0:
                codon2_frame1.append(index)
            elif index % 3 == 1:
                codon2_frame2.append(index)
            elif index % 3 == 2:
                codon2_frame3.append(index)
            else:
                pass

        # if any of the codon2 indexes are higher than those in the codon1 for each frame, or if codon one frame is empty, change the frames to complete
        if any(a > b for a in codon2_frame1 for b in codon1_frame1) or not codon1_frame1:
            self.frame1_complete = True
        if any(a > b for a in codon2_frame2 for b in codon1_frame2) or not codon1_frame2:
            self.frame2_complete = True
        if any(a > b for a in codon2_frame3 for b in codon1_frame3) or not codon1_frame3:
            self.frame3_complete = True
        #if all frames are complete, set all_frames_complete to true
        if self.frame1_complete == True and self.frame2_complete == True and self.frame3_complete == True:
            self.all_frames_complete = True

    #merge paths between nodes, and update the orientiation of the final node
    def merge_path(self, self_seq, self_dir, new_node, new_seq_dir, ksize, GFA):
        merged_path = None
        new_seq = GFA._graph.node[str(new_node)]['sequence']
        if self_dir == "+" and new_seq_dir == "+":
            merged_path = self_seq + new_seq[(ksize - 1)::]
            self.ori = "+"
        elif self_dir == "+" and new_seq_dir == "-":
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = self_seq + rc_new_seq[(ksize - 1)::]
            self.ori = "-"
        elif self_dir == "-" and new_seq_dir == "+":
            seq1 = Seq(self.seq)
            rc_seq1 = str(seq1.reverse_complement())
            merged_path = new_seq + rc_seq1[(ksize - 1)::]
            self.ori = "+"
        elif self_dir == "-" and new_seq_dir == "-":
            seq1 = Seq(self_seq)
            rc_seq1 = str(seq1.reverse_complement())
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = rc_new_seq + rc_seq1[(ksize - 1)::]
            self.ori = "-"
        else:
            pass
        return merged_path

class Path_old:
    def __init__(self, GFA, nodes, ksize, codon1=None, codon2=None, frame1_complete=False, frame2_complete=False, frame3_complete=False):
        #initialising entries for path, where nodes is a list
        self.nodes = [nodes[0]]

        #initialise orientation; absori being orientation of first node in whole path, relori being orientation of last node added
        self.absori = "+"
        self.relori = "+"

        #create edge list for beginning node
        self.edges = self.update_edges(self.nodes[0], GFA)

        #add sequences for detailed nodes if they are connected to eachother. If repeat = False, ignore target nodes already in path and set all frame complete to true
        self.seq = GFA._graph.node[str(self.nodes[0])]['sequence']
        if len(nodes) == 1:
            pass
        else:
            for index, i in enumerate(nodes[1:]):
                edge_list = self.edges
                for edge in edge_list:
                    source, source_dir, target, target_dir = edge
                    #merge sequences isn't working, need to take into account the orientiation of the whole merged unitig!
                    if source == str(nodes[index]) and target == str(i):
                        self.seq = self.merge_path(self.seq, source_dir, target, target_dir, ksize, GFA)
                        self.nodes.append(target)
                        self.edges = self.update_edges(self.nodes[-1], GFA)
                        self.absori = source_dir
                        self.relori = target_dir
                    else:
                        pass
        self.len = len(self.seq)

        #generate codon k-mers for path
        kmer_list = []
        for i in range(0, (len(self.seq) - 2)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

        #create booleans to check if all paths are complete for iterative algorithm, update if codons are present
        self.frame1_complete = frame1_complete
        self.frame2_complete = frame2_complete
        self.frame3_complete = frame3_complete

        #check if k-mers are present depending on absolute orientiation of merged unitig based on first unitig, reversing ordering of codon indexing.
        if frame1_complete == False:
            if self.absori == "+":
                self.check_frame(codon1, codon2, 1, False)
            elif self.absori == "-":
                rc_codon1 = str(Seq(codon1).reverse_complement())
                rc_codon2 = str(Seq(codon2).reverse_complement())
                self.check_frame(rc_codon1, rc_codon2, 1, True)

        if frame2_complete == False:
            if self.absori == "+":
                self.check_frame(codon1, codon2, 2, False)
            elif self.absori == "-":
                rc_codon1 = str(Seq(codon1).reverse_complement())
                rc_codon2 = str(Seq(codon2).reverse_complement())
                self.check_frame(rc_codon1, rc_codon2, 2, True)

        if frame3_complete == False:
            if self.absori == "+":
                self.check_frame(codon1, codon2, 3, False)
            elif self.absori == "-":
                rc_codon1 = str(Seq(codon1).reverse_complement())
                rc_codon2 = str(Seq(codon2).reverse_complement())
                self.check_frame(rc_codon1, rc_codon2, 3, True)

        #check if all paths are complete
        self.all_frames_complete = False
        self.check_all_complete()

    #class method to genereate edges from input node
    def update_edges(self, node, GFA):
        edge_list = []
        for sink_nodeid, sink_node_info in GFA._graph.edge[str(node)].items():
            # iterate through node info, picking out one of either virtual edges where the node of interest is the 'from_orn'
            for virtual_edgeid, virtual_edge_info in sink_node_info.items():
                if virtual_edge_info['from_node'] == str(node):
                    node_tuple = (
                        virtual_edge_info['from_node'], virtual_edge_info['from_orn'], virtual_edge_info['to_node'],
                        virtual_edge_info['to_orn'])
                    edge_list.append(node_tuple)
        return edge_list

    # class method to enable updating of kmers in seq from last node merged
    def updates_kmers(self):
        kmer_list = []
        for i in range(0, (len(self.seq) - 3)):
            kmer = self.seq[i:i + 3]
            kmer_list.append(kmer)
        self.codons = kmer_list

    #check designated frame to see if it is complete. Also check if codon1 and codon2 are reverse due to reverse complementation
    def check_frame(self, codon1, codon2, frame, reversed):
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
        if reversed == False:
            if any(a > b for a in codon2_frame for b in codon1_frame) or not codon1_frame:
                if frame == 1:
                    self.frame1_complete = True
                elif frame == 2:
                    self.frame2_complete = True
                elif frame == 3:
                    self.frame3_complete = True
                else:
                    pass
        elif reversed == True:
            if any(a < b for a in codon2_frame for b in codon1_frame) or not codon1_frame:
                if frame == 1:
                    self.frame1_complete = True
                elif frame == 2:
                    self.frame2_complete = True
                elif frame == 3:
                    self.frame3_complete = True
                else:
                    pass
        else:
            pass

    def check_all_complete(self):
        if self.frame1_complete == True and self.frame2_complete == True and self.frame3_complete == True:
            self.all_frames_complete = True

    #merge paths between nodes, and update the orientiation of the final node
    def merge_path(self, self_seq, self_dir, new_node, new_seq_dir, ksize, GFA):
        merged_path = None
        new_seq = GFA._graph.node[str(new_node)]['sequence']
        if self_dir == "+" and new_seq_dir == "+":
            merged_path = self_seq + new_seq[(ksize - 1)::]
        elif self_dir == "+" and new_seq_dir == "-":
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = self_seq + rc_new_seq[(ksize - 1)::]
        elif self_dir == "-" and new_seq_dir == "+":
            seq1 = Seq(self_seq)
            rc_seq1 = str(seq1.reverse_complement())
            merged_path = rc_seq1 + new_seq[(ksize - 1)::]
        elif self_dir == "-" and new_seq_dir == "-":
            seq1 = Seq(self_seq)
            rc_seq1 = str(seq1.reverse_complement())
            new_seq = Seq(new_seq)
            rc_new_seq = str(new_seq.reverse_complement())
            merged_path = rc_seq1 + rc_new_seq[(ksize - 1)::]
        else:
            pass
        return merged_path

# this method aims to merge unitigs and count if there are stop codons in each of the three reading frames, then just count whether there are one of each codon in correct orientation
def iter_codons_merge(GFA, nodes_list, codon1, codon2, ksize):
    # create all gene paths to add paths to
    all_gene_paths = {}

    for i in nodes_list:
        #create paths list to hold all complete and incomplete paths
        complete_path_list = []
        incomplete_path_list = []

        #create intial path, id 0
        init_path = Path(0, i, GFA)

        # index k-mers based on chosen codons
        init_path.index_kmers(codon1, codon2)

        #check to see if initial unitig has complete paths
        if init_path.all_frames_complete == True:
            complete_path_list.append(init_path)
            all_paths_complete = True
        else:
            incomplete_path_list.append(init_path)
            all_paths_complete = False

        #excute while loop to iteratively go through all paths generated by edges
        path_id = 1
        while all_paths_complete == False:
            iter_path_list = []
            for path in incomplete_path_list:
                #generate a new path for each edge
                for edge in path.edges:
                    #unpack path end edges
                    source, source_dir, target, target_dir = edge

                    #prevent same path revisiting nodes already visited
                    if target not in path.nodes:
                        # create full copy of original path, reset frame complete status
                        new_path = copy.deepcopy(path)
                        new_path.change_id(path_id)

                        path_id += 1
                        new_path.merge_path(source_dir, target, target_dir, ksize, GFA)
                        new_path.reset_frame()
                        new_path.index_kmers(codon1, codon2)

                        ### for debugging ###
                        new_path_seq = new_path.seq
                        new_path_codons = new_path.codons
                        new_path_frame1 = new_path.frame1_complete
                        new_path_frame2 = new_path.frame2_complete
                        new_path_frame3 = new_path.frame3_complete
                        new_path_nodes = new_path.nodes
                        new_path_all_frames = new_path.all_frames_complete

                        # check if path has complete frames, else add it to the iter_path_list
                        if new_path.all_frames_complete == True:
                            complete_path_list.append(new_path)
                        else:
                            iter_path_list.append(new_path)
                    else:
                        pass

                #check if iter_path_list is empty, this means that no paths are incomplete. If not empty, run index_kmers class method
                #also empty incomplete path list, as all paths are either stored in complete_path list or iter_path_list
                incomplete_path_list = []
                if not iter_path_list:
                    all_paths_complete = True
                else:
                    for path in iter_path_list:
                        path.index_kmers(codon1, codon2)
                        incomplete_path_list.append(path)

        #return complete paths for selected node
        all_gene_paths[i] = [str(item.seq) for item in complete_path_list]
    return all_gene_paths

#recursive algorithm to generate strings and nodes with codons
def recur_paths(GFA, start_node_list, codon1, codon2, ksize):

    path_list = [start_node_list]
    start_path = Path2(GFA, start_node_list[-1], start_node_list, ksize, codon1=codon1, codon2=codon2)

    if start_path.all_frames_complete == False:
        for neighbour in start_path.edges:
            source, source_dir, target, target_dir = neighbour
            path = start_node_list + [target]
            for iteration in recur_paths(GFA, path, codon1, codon2, ksize):
                path_list.append(iteration)
        del path_list[0]
    else:
        pass
    return path_list


###for debugging###

if __name__ == '__main__':
    from pygfa import *
    import networkx as nx
    from Bio.Seq import Seq
    import copy
    graph = pygfa.gfa.GFA.from_file("test.gfa")

    node_list = [1, 2, 3]
    path1 = Path2(graph, "1", node_list, 31)
    #iter_codons_merge(graph, nodes_stop, "ACT", "ACT", 31)