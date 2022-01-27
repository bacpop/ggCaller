from collections import defaultdict, Counter
from Bio.Seq import translate, reverse_complement, Seq
from multiprocessing import get_context
from Bio import SeqIO
from panaroo.cdhit import align_dna_cdhit
from panaroo.isvalid import del_dups
from joblib import Parallel, delayed
import os
import gffutils as gff
from io import StringIO
import edlib
from panaroo.merge_nodes import delete_node, remove_member_from_node
from tqdm import tqdm
import re

from ggCaller.shared_memory import *
from functools import partial


def find_missing(G,
                 graph_shd_arr_tup,
                 high_scoring_ORFs,
                 is_ref,
                 write_idx,
                 kmer,
                 repeat,
                 isolate_names,
                 merge_id_thresh,
                 search_radius,
                 prop_match,
                 pairwise_id_thresh,
                 n_cpu,
                 remove_by_consensus=False,
                 verbose=True):
    # load shared memory items
    graph_existing_shm = shared_memory.SharedMemory(name=graph_shd_arr_tup.name)
    graph_shd_arr = np.ndarray(graph_shd_arr_tup.shape, dtype=graph_shd_arr_tup.dtype, buffer=graph_existing_shm.buf)

    # iterate through nodes to identify accessory genes for searching
    # these are nodes missing a member with at least one neighbour that has that member
    n_searches = 0
    search_dict = {key: {"conflicts": {}, "searches": {}} for key in range(0, len(isolate_names))}
    for node in G.nodes():
        for neigh in G.neighbors(node):
            for sid in sorted(G.nodes[neigh]['seqIDs']):
                member = int(sid.split("_")[0])
                ORF_ID = int(sid.split("_")[-1])
                ORF_info = high_scoring_ORFs[member][ORF_ID]
                search_dict[member]["conflicts"][neigh] = ORF_info

                if member not in G.nodes[node]['members']:
                    if len(G.nodes[node]["dna"][G.nodes[node]
                    ['maxLenId']]) <= 0:
                        print(G.nodes[node]["dna"])
                        raise NameError("Problem!")
                    # add the representative DNA sequence for missing node and the ID of the colour to search from
                    if node not in search_dict[member]["searches"]:
                        search_dict[member]["searches"][node] = {}
                    search_dict[member]["searches"][node][G.nodes[node]["dna"][G.nodes[node]['maxLenId']]] = ORF_info[
                                                                                                             :6]

                    n_searches += 1

    if verbose:
        print("Number of searches to perform: ", n_searches)
        print("Searching...")

    # create lists to return results
    all_hits = [None] * len(isolate_names)
    all_node_locs = [None] * len(isolate_names)

    with Pool(processes=n_cpu, maxtasksperchild=1) as pool:
        for member, hits, node_locs in pool.map(partial(search_graph,
                                                        graph_shd_arr_tup=graph_shd_arr_tup,
                                                        isolate_names=isolate_names,
                                                        search_radius=search_radius,
                                                        prop_match=prop_match,
                                                        pairwise_id_thresh=pairwise_id_thresh,
                                                        is_ref=is_ref,
                                                        kmer=kmer,
                                                        repeat=repeat),
                                                search_dict.items()):
            all_hits[member] = hits
            all_node_locs[member] = node_locs

    if verbose:
        print("translating hits...")

    hits_trans_dict = {}
    for member, hits in enumerate(all_hits):
        hits_trans_dict[member] = Parallel(n_jobs=n_cpu)(
            delayed(translate_to_match)(hit[1], G.nodes[hit[0]]["protein"][0])
            for hit in hits)

    # remove nodes that conflict (overlap)
    nodes_by_size = sorted([(G.nodes[node]['size'], node)
                            for node in G.nodes()],
                           reverse=True)
    nodes_by_size = [n[1] for n in nodes_by_size]
    bad_node_mem_pairs = set()
    bad_nodes = set()
    for member, node_locs in enumerate(all_node_locs):
        seq_coverage = {}
        # iterate over all nodes and add sequence information
        for node in nodes_by_size:
            # if node in bad_nodes: continue
            if node not in node_locs: continue
            # iterate over all the nodes present for the current ORF
            node_coverage = 0
            for DBG_node, node_coords in zip(node_locs[node][0], node_locs[node][1]):
                # check if node is negative and needs reversing
                # make copies to avoid editing in place
                temp_DBG_node = DBG_node
                if temp_DBG_node < 0:
                    temp_DBG_node *= -1
                    node_end = graph_shd_arr[0].node_size(temp_DBG_node) - 1
                    temp_node_coords = (node_end - node_coords[1], node_end - node_coords[0])
                else:
                    temp_node_coords = node_coords
                if temp_DBG_node not in seq_coverage:
                    seq_coverage[temp_DBG_node] = np.zeros(temp_node_coords[1] + 1, dtype=bool)
                else:
                    node_coverage += np.sum(seq_coverage[temp_DBG_node][temp_node_coords[0]:temp_node_coords[1]])
            # negate coverage from node overlaps in DBG
            node_coverage -= node_locs[node][2]
            # check if node coverage is exceeded
            if node_coverage >= 0.5 * (max(G.nodes[node]['lengths'])):
                if member in G.nodes[node]['members']:
                    remove_member_from_node(G, node, member)
                bad_node_mem_pairs.add((node, member))
            else:
                # if sequence coverage is less than threshold, iterate again and add sequence information
                for DBG_node, node_coords in zip(node_locs[node][0], node_locs[node][1]):
                    temp_DBG_node = DBG_node
                    if temp_DBG_node < 0:
                        temp_DBG_node *= -1
                        node_end = graph_shd_arr[0].node_size(temp_DBG_node) - 1
                        temp_node_coords = (node_end - node_coords[1], node_end - node_coords[0])
                    else:
                        temp_node_coords = node_coords
                    seq_coverage[temp_DBG_node][temp_node_coords[0]:temp_node_coords[1]] = True

    for node in G.nodes():
        if len(G.nodes[node]['members']) <= 0:
            bad_nodes.add(node)
    for node in bad_nodes:
        if node in G.nodes():
            delete_node(G, node)

    # remove by consensus
    if remove_by_consensus:
        if verbose:
            print("removing by consensus...")
        node_hit_counter = Counter()
        for member, hits in enumerate(all_hits):
            for node, dna_hit in hits:
                if dna_hit == "": continue
                if node in bad_nodes: continue
                if (node, member) in bad_node_mem_pairs: continue
                node_hit_counter[node] += 1
        for node in G:
            if node_hit_counter[node] > G.nodes[node]['size']:
                bad_nodes.add(node)
        for node in bad_nodes:
            if node in G.nodes():
                delete_node(G, node)

    if verbose:
        print("Updating output...")

    # update the graph
    n_found = 0
    for member, hits in enumerate(all_hits):
        i = -1
        for node, dna_hit in hits:
            i += 1
            if dna_hit == "": continue
            if node in bad_nodes: continue
            if (node, member) in bad_node_mem_pairs: continue
            n_found += 1
            hit_protein, hit_dna = hits_trans_dict[member][i]
            G.nodes[node]['members'].add(member)
            G.nodes[node]['size'] += 1
            G.nodes[node]['dna'].append(dna_hit)
            G.nodes[node]['protein'].append(hit_protein)
            G.nodes[node]['seqIDs'] |= set([str(member) + "_refound_" + str(n_found * -1)])
            # add new refound gene to high_scoring_ORFs with negative ID to indicate refound
            nodelist, node_coords, total_overlap = all_node_locs[member][node]
            premature_stop = "*" in hit_protein[1:-3]
            if G.nodes[node]['bitscore'] != 0:
                annotation = ("refound", G.nodes[node]['annotation'],
                              G.nodes[node]['bitscore'], G.nodes[node]['description'])
            else:
                annotation = ("refound", "hypothetical protein", 0, "hypothetical protein")
            high_scoring_ORFs[member][n_found * -1] = (
                nodelist, node_coords, len(dna_hit), premature_stop, hit_protein, hit_dna, annotation)

    if verbose:
        print("Number of refound genes: ", n_found)

    return G, high_scoring_ORFs


def search_graph(search_pair,
                 graph_shd_arr_tup,
                 isolate_names,
                 is_ref,
                 kmer,
                 repeat,
                 search_radius=10000,
                 prop_match=0.2,
                 pairwise_id_thresh=0.95):
    # unpack search_pair and assign fasta
    member, dicts = search_pair
    fasta = isolate_names[member]

    # load shared memory items
    graph_existing_shm = shared_memory.SharedMemory(name=graph_shd_arr_tup.name)
    graph_shd_arr = np.ndarray(graph_shd_arr_tup.shape, dtype=graph_shd_arr_tup.dtype, buffer=graph_existing_shm.buf)

    # sort items to preserve order
    conflicts = {k: v for k, v in sorted(dicts["conflicts"].items(), key=lambda item: item[0])}
    node_search_dict = dicts["searches"]

    node_locs = {}

    # mask regions that already have genes
    for node, ORF_info in conflicts.items():
        # determine sequence overlap of ORFs
        total_overlap = 0
        for i, node_coords in enumerate(ORF_info[1]):
            if i != 0:
                if node_coords[1] >= kmer - 1:
                    total_overlap += ((kmer - 1) - node_coords[0])
                else:
                    total_overlap += (node_coords[1] - node_coords[0]) + 1
        node_locs[node] = (ORF_info[0], ORF_info[1], total_overlap)

    # get sequences to search
    node_search_dict = graph_shd_arr[0].refind_gene(member, node_search_dict, search_radius, is_ref,
                                                    kmer, fasta, repeat)

    # search for matches
    hits = []
    for node in node_search_dict:
        best_hit = ""
        best_loc = None
        for search, search_info in node_search_dict[node].items():
            db_seq, nodelist, node_ranges, path_rev_comp = search_info

            # check if no sequence found
            if db_seq == "":
                pass

            hit, loc, rev_comp = search_dna(db_seq,
                                            search,
                                            prop_match,
                                            pairwise_id_thresh,
                                            refind=True)

            # work out how ORF is orientated with respect to original sequence (only if FM-index used)
            if rev_comp and path_rev_comp and is_ref:
                rev_comp = False

            # convert linear coordinates into node coordinates
            ORF_graph_loc = convert_coords(loc, nodelist, node_ranges, kmer - 1)

            # if db_seq was reversed to align, need to reverse node coordinates
            if rev_comp and hit != "":
                reversed_nodes = []
                reversed_loci = []
                for node_ID, node_coords in zip(ORF_graph_loc[0], ORF_graph_loc[1]):
                    rev_node_ID = node_ID * -1
                    node_end = graph_shd_arr[0].node_size(rev_node_ID) - 1
                    reversed_end = node_end - node_coords[0]
                    reversed_start = node_end - node_coords[1]
                    reversed_nodes.append(rev_node_ID)
                    reversed_loci.append((reversed_start, reversed_end))
                reversed_nodes.reverse()
                reversed_loci.reverse()
                ORF_graph_loc = (reversed_nodes, reversed_loci, ORF_graph_loc[-1])


            if len(hit) > len(best_hit):
                best_hit = hit
                best_loc = ORF_graph_loc

        hits.append((node, best_hit))
        if (best_loc is not None) and (best_hit != ""):
            node_locs[node] = best_loc

    return member, hits, node_locs


def convert_coords(loc, nodelist, node_ranges, overlap):
    # iterate over nodelist and node_coords to determine what bases are traversed
    traversed_nodes = []
    traversed_loci = []

    total_overlap = 0

    # check if hit found, otherwise pass
    if loc[1] != 0:
        # adjust loc[1] due to different indexing in alignment
        loc[1] = loc[1] - 1
        for i in range(0, len(nodelist)):
            start_assigned = False
            end_assigned = False

            # if start of ORF is below node range, then check ORF traverses node
            if loc[0] < node_ranges[i][0]:
                traversed_node_start = 0
                start_assigned = True
            elif loc[0] >= node_ranges[i][0] and loc[0] < node_ranges[i][1]:
                traversed_node_start = loc[0] - node_ranges[i][0]
                start_assigned = True

            # if end of ORF is above node range, then check if ORF traversed
            if loc[1] >= node_ranges[i][1]:
                traversed_node_end = node_ranges[i][2]
                end_assigned = True
            elif loc[1] >= node_ranges[i][0] and loc[1] < node_ranges[i][1]:
                traversed_node_end = loc[1] - node_ranges[i][0]
                end_assigned = True

            # if the ORF traverses node, update coordinates
            if start_assigned and end_assigned:
                traversed_nodes.append(nodelist[i])
                traversed_loci.append((traversed_node_start, traversed_node_end))

                # work out how much sequence overlap is present between nodes, add 1 and zero indexed
                if len(traversed_nodes) != 1:
                    if traversed_node_end >= overlap:
                        total_overlap += (overlap - traversed_node_start)
                    else:
                        total_overlap += (traversed_node_end - traversed_node_start) + 1


            # gone past last node covering ORF so assign 3p and end index to previous node
            elif start_assigned and not end_assigned:
                break

    return traversed_nodes, traversed_loci, total_overlap


def repl(m):
    return ('X' * len(m.group()))


def search_dna(db_seq, search_sequence, prop_match, pairwise_id_thresh,
               refind):
    found_dna = ""
    start = None
    end = None
    max_hit = 0
    loc = [0, 0]
    rev_comp = False

    added_E_len = int(len(search_sequence) / 2)

    for i, db in enumerate([db_seq, str(Seq(db_seq).reverse_complement())]):

        # add some Ns at the start and end to deal with fragments at the end of contigs
        db = "E" * added_E_len + db + "E" * added_E_len

        aln = edlib.align(search_sequence,
                          db,
                          mode="HW",
                          task='path',
                          k=10 * len(search_sequence),
                          additionalEqualities=[
                              ('A', 'N'),
                              ('C', 'N'),
                              ('G', 'N'),
                              ('T', 'N'),
                              ('A', 'E'),
                              ('C', 'E'),
                              ('G', 'E'),
                              ('T', 'E'),
                          ])

        # remove trailing inserts
        cig = re.split(r'(\d+)', aln['cigar'])[1:]
        if cig[-1] == "I":
            aln['editDistance'] -= int(cig[-2])
        if cig[1] == "I":
            aln['editDistance'] -= int(cig[0])

        if aln['editDistance'] == -1:
            start = -1
        else:
            # take hit that is closest to the centre of the neighbouring gene
            centre = len(db) / 2.0
            tloc = min(aln['locations'],
                       key=lambda x: min(centre - x[0], centre - x[1]))
            start = tloc[0]
            end = tloc[1] + 1

        # if found:
        #     print(aln)
        #     print(start, end)

        # skip if nothing was found
        if start == -1: continue

        possible_dbs = [db]
        if db.find("NNNNNNNNNNNNNNNNNNNN") != -1:
            possible_dbs += [
                re.sub("^[ACGTEX]{0,}NNNNNNNNNNNNNNNNNNNN", repl, db, 1),
                re.sub("NNNNNNNNNNNNNNNNNNNN[ACGTEX]{0,}$", repl, db, 1)
            ]

        for posdb in possible_dbs:
            # skip if alignment is too short
            n_X = posdb[start:end].count("X")
            n_E = posdb[start:end].count("E")

            aln_length = float(end - start - n_X - n_E)
            if (aln_length / len(search_sequence)) <= prop_match: continue
            if (posdb[start:end].count("A") + posdb[start:end].count("C") +
                posdb[start:end].count("G") + posdb[start:end].count("T")
            ) / len(search_sequence) <= prop_match:
                continue

            # determine an approximate percentage identity
            pid = 1.0 - (aln['editDistance'] - n_X) / (1.0 * aln_length)

            # skip if identity below threshold
            if pid <= pairwise_id_thresh: continue

            # if found:
            #     print("aln_length:", aln_length)
            #     print("pid:", pid)

            if max_hit < (pid * aln_length):
                found_dna = posdb[start:end]
                max_hit = (pid * aln_length)
                if i == 0:
                    loc = [start, end]
                else:
                    loc = [len(posdb) - tloc[1] - 1, len(posdb) - tloc[0]]
                    rev_comp = True
                loc = [
                    max(0,
                        min(loc) - added_E_len),
                    min(max(loc) - added_E_len, len(db_seq))
                ]

    seq = found_dna.replace('X', 'N').replace('E', 'N')
    seq = seq.strip('N')

    return seq, loc, rev_comp


def translate_to_match(hit, target_prot):
    if hit == "": return ""

    # translate in all 6 frames splitting on unknown
    dna_seqs = [
        s[i:].ljust(len(s[i:]) + (3 - len(s[i:]) % 3), 'N')
        for i in range(3) for s in [hit, reverse_complement(hit)]
    ]

    proteins = [
        translate(s) for s in dna_seqs
    ]

    search_set = set(
        [target_prot[i:i + 3] for i in range(len(target_prot) - 2)])

    alignments = []
    for index, target_sequence in enumerate(proteins):
        query_set = set([
            target_sequence[i:i + 3] for i in range(len(target_sequence) - 2)
        ])
        alignments.append(
            (target_sequence, len(search_set.intersection(query_set)), index))

    prot = max(alignments, key=lambda x: x[1])

    return prot[0], dna_seqs[prot[2]]


blosum50 = \
    {
        '*': {'*': 1, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5,
              'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5,
              'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5,
              'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5},
        'A': {'*': -5, 'A': 5, 'C': -1, 'B': -2, 'E': -1, 'D': -2, 'G': 0,
              'F': -3, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -2,
              'N': -1, 'Q': -1, 'P': -1, 'S': 1, 'R': -2, 'T': 0, 'W': -3,
              'V': 0, 'Y': -2, 'X': -1, 'Z': -1},
        'C': {'*': -5, 'A': -1, 'C': 13, 'B': -3, 'E': -3, 'D': -4,
              'G': -3, 'F': -2, 'I': -2, 'H': -3, 'K': -3, 'M': -2,
              'L': -2, 'N': -2, 'Q': -3, 'P': -4, 'S': -1, 'R': -4,
              'T': -1, 'W': -5, 'V': -1, 'Y': -3, 'X': -1, 'Z': -3},
        'B': {'*': -5, 'A': -2, 'C': -3, 'B': 6, 'E': 1, 'D': 6, 'G': -1,
              'F': -4, 'I': -4, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 5,
              'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -5, 'V': -3,
              'Y': -3, 'X': -1, 'Z': 1},
        'E': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 6, 'D': 2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': -1, 'R': 0, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5},
        'D': {'*': -5, 'A': -2, 'C': -4, 'B': 6, 'E': 2, 'D': 8, 'G': -1,
              'F': -5, 'I': -4, 'H': -1, 'K': -1, 'M': -4, 'L': -4, 'N': 2,
              'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -5, 'V': -4,
              'Y': -3, 'X': -1, 'Z': 1},
        'G': {'*': -5, 'A': 0, 'C': -3, 'B': -1, 'E': -3, 'D': -1, 'G': 8,
              'F': -4, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0,
              'Q': -2, 'P': -2, 'S': 0, 'R': -3, 'T': -2, 'W': -3, 'V': -4,
              'Y': -3, 'X': -1, 'Z': -2},
        'F': {'*': -5, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -5,
              'G': -4, 'F': 8, 'I': 0, 'H': -1, 'K': -4, 'M': 0, 'L': 1,
              'N': -4, 'Q': -4, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 1,
              'V': -1, 'Y': 4, 'X': -1, 'Z': -4},
        'I': {'*': -5, 'A': -1, 'C': -2, 'B': -4, 'E': -4, 'D': -4,
              'G': -4, 'F': 0, 'I': 5, 'H': -4, 'K': -3, 'M': 2, 'L': 2,
              'N': -3, 'Q': -3, 'P': -3, 'S': -3, 'R': -4, 'T': -1,
              'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -3},
        'H': {'*': -5, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2,
              'F': -1, 'I': -4, 'H': 10, 'K': 0, 'M': -1, 'L': -3, 'N': 1,
              'Q': 1, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -3, 'V': -4,
              'Y': 2, 'X': -1, 'Z': 0},
        'K': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 6, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': 0, 'R': 3, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 1},
        'M': {'*': -5, 'A': -1, 'C': -2, 'B': -3, 'E': -2, 'D': -4,
              'G': -3, 'F': 0, 'I': 2, 'H': -1, 'K': -2, 'M': 7, 'L': 3,
              'N': -2, 'Q': 0, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -1,
              'V': 1, 'Y': 0, 'X': -1, 'Z': -1},
        'L': {'*': -5, 'A': -2, 'C': -2, 'B': -4, 'E': -3, 'D': -4,
              'G': -4, 'F': 1, 'I': 2, 'H': -3, 'K': -3, 'M': 3, 'L': 5,
              'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -1,
              'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3},
        'N': {'*': -5, 'A': -1, 'C': -2, 'B': 5, 'E': 0, 'D': 2, 'G': 0,
              'F': -4, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -4, 'N': 7,
              'Q': 0, 'P': -2, 'S': 1, 'R': -1, 'T': 0, 'W': -4, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 0},
        'Q': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2,
              'F': -4, 'I': -3, 'H': 1, 'K': 2, 'M': 0, 'L': -2, 'N': 0,
              'Q': 7, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -1, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 4},
        'P': {'*': -5, 'A': -1, 'C': -4, 'B': -2, 'E': -1, 'D': -1,
              'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -3,
              'L': -4, 'N': -2, 'Q': -1, 'P': 10, 'S': -1, 'R': -3,
              'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': -1},
        'S': {'*': -5, 'A': 1, 'C': -1, 'B': 0, 'E': -1, 'D': 0, 'G': 0,
              'F': -3, 'I': -3, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1,
              'Q': 0, 'P': -1, 'S': 5, 'R': -1, 'T': 2, 'W': -4, 'V': -2,
              'Y': -2, 'X': -1, 'Z': 0},
        'R': {'*': -5, 'A': -2, 'C': -4, 'B': -1, 'E': 0, 'D': -2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 3, 'M': -2, 'L': -3, 'N': -1,
              'Q': 1, 'P': -3, 'S': -1, 'R': 7, 'T': -1, 'W': -3, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 0},
        'T': {'*': -5, 'A': 0, 'C': -1, 'B': 0, 'E': -1, 'D': -1, 'G': -2,
              'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0,
              'Q': -1, 'P': -1, 'S': 2, 'R': -1, 'T': 5, 'W': -3, 'V': 0,
              'Y': -2, 'X': -1, 'Z': -1},
        'W': {'*': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -3, 'D': -5,
              'G': -3, 'F': 1, 'I': -3, 'H': -3, 'K': -3, 'M': -1, 'L': -2,
              'N': -4, 'Q': -1, 'P': -4, 'S': -4, 'R': -3, 'T': -3,
              'W': 15, 'V': -3, 'Y': 2, 'X': -1, 'Z': -2},
        'V': {'*': -5, 'A': 0, 'C': -1, 'B': -3, 'E': -3, 'D': -4, 'G': -4,
              'F': -1, 'I': 4, 'H': -4, 'K': -3, 'M': 1, 'L': 1, 'N': -3,
              'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 5,
              'Y': -1, 'X': -1, 'Z': -3},
        'Y': {'*': -5, 'A': -2, 'C': -3, 'B': -3, 'E': -2, 'D': -3,
              'G': -3, 'F': 4, 'I': -1, 'H': 2, 'K': -2, 'M': 0, 'L': -1,
              'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -1, 'T': -2, 'W': 2,
              'V': -1, 'Y': 8, 'X': -1, 'Z': -2},
        'X': {'*': -5, 'A': -1, 'C': -1, 'B': -1, 'E': -1, 'D': -1,
              'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1,
              'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1,
              'T': -1, 'W': -1, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1},
        'Z': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 5, 'D': 1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0,
              'Q': 4, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -2, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5}}
