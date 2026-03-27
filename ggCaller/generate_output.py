from functools import partial
from ggCaller.shared_memory import *
import ggCaller_cpp
from ggCaller import __version__
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import gzip
from functools import partial

# stringizer imports
import os
from ast import literal_eval
from collections import defaultdict

try:
    long
except NameError:
    long = int
try:
    unicode
except NameError:
    unicode = str
try:
    unichr
except NameError:
    unichr = chr
try:
    literal_eval(r"u'\u4444'")
except SyntaxError:
    # Remove 'u' prefixes in unicode literals in Python 3
    def rtp_fix_unicode(s):
        return s[1:]
else:
    rtp_fix_unicode = None


def print_ORF_calls(ORF_file_paths, outfile, input_colours, overlap, DBG, truncation_threshold=0, G=None, ids_len_stop=None):
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]
    DNA_records = []
    aa_records = []
    if G == None or ids_len_stop == None:
        for colour, file_path in ORF_file_paths.items():
            ORF_map = ggCaller_cpp.read_ORF_file(file_path)

            for ORF_ID, ORFNodeVector in ORF_map.items():
                gene = DBG.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
                DNA_records.append(SeqRecord(Seq(gene), id=isolate_names[colour] + "_" + str(ORF_ID),
                                             description=""))
                aa_records.append(SeqRecord(Seq(gene).translate(), id=isolate_names[colour] + "_" + str(ORF_ID),
                                            description=""))
    else:
        for colour, file_path in ORF_file_paths.items():
            to_remove = set()
            ORF_map = ggCaller_cpp.read_ORF_file(file_path)

            for ORF_ID, ORF_info in ORF_map.items():
                delim = "_0_" if ORF_ID >= 0 else "_refound_"
                pan_ORF_id = str(colour) + delim + str(ORF_ID)
                # if pan_ORF_id not in ORF_map means has been removed so remove it from ORF_map
                node = ids_len_stop[pan_ORF_id][2]
                if node == -1:
                    to_remove.add(ORF_ID)
                    del ids_len_stop[pan_ORF_id]
                    continue
                node_annotation = G.nodes[node]["description"]
                length_centroid = G.nodes[node]['lengths'][G.nodes[node]['maxLenId']]

                # generate output sequences
                gene = DBG.generate_sequence(ORF_info[0], ORF_info[1], overlap)
                ORF_len = ORF_info[2]
                gene_annotation = node_annotation
                if ORF_len < (length_centroid * truncation_threshold) or (
                        ORF_ID < 0 and (ORF_info[3] is True
                                        or ORF_len % 3 != 0)):
                    gene_annotation += " (potential psuedogene)"
                DNA_records.append(SeqRecord(Seq(gene), id=isolate_names[colour] + "_" + str(ORF_ID),
                                description="ggcID=" + pan_ORF_id + ";" + gene_annotation))
                if (ORF_ID < 0):
                    aa_records.append(
                        SeqRecord(Seq(ORF_info[5]), id=isolate_names[colour] + "_" + str(ORF_ID),
                                        description="ggcID=" + pan_ORF_id + ";" + gene_annotation))
                else:
                    aa_records.append(
                        SeqRecord(Seq(gene).translate(), id=isolate_names[colour] + "_" + str(ORF_ID),
                                        description="ggcID=" + pan_ORF_id + ";" + gene_annotation))
            
            # remove from ORF_map and save
            for ORF_ID in to_remove:
                del ORF_map[ORF_ID]
            ggCaller_cpp.save_ORF_file(file_path, ORF_map)
    SeqIO.write(DNA_records, outfile + ".ffn", "fasta")
    SeqIO.write(aa_records, outfile + ".faa", "fasta")
    return

def generate_GFF(graph, ORF_file_paths, input_colours, isolate_names, contig_annotation, output_dir, overlap,
                 write_idx, ref_list, num_threads, Path_dir):
    # create directory for gffs
    GFF_dir = os.path.join(output_dir, "GFF")
    if not os.path.exists(GFF_dir):
        os.mkdir(GFF_dir)
    # make sure trailing forward slash is present
    GFF_dir = os.path.join(GFF_dir, "")

    # iterate over colours in contig_annotation, writing output file
    for colour, file_path in ORF_file_paths.items():
        # if colour is not reference, then cannot generate GFF
        if not ref_list[colour]:
            continue
        
        ORF_map = ggCaller_cpp.read_ORF_file(file_path)

        # get coordinates for each ORF
        ORF_IDs = [(ORF_map[i[0]][0], ORF_map[i[0]][1]) for i in
                   contig_annotation[colour]]
        ORF_coords = graph.ORF_location(ORF_IDs, input_colours[colour], overlap, write_idx, num_threads, Path_dir)

        # create data structure for each GFF entry
        GFF_entries = defaultdict(list)
        for i in range(len(contig_annotation[colour])):
            GFF_entries[ORF_coords[i][0][0]].append((contig_annotation[colour][i][0],
                                                     (ORF_coords[i][0][1], ORF_coords[i][1]),
                                                     contig_annotation[colour][i][1]))

        # check if FASTA is gzipped
        gzipped = False
        with open(input_colours[colour], 'rb') as test_f:
            gzipped = True if test_f.read(2) == b'\x1f\x8b' else False

        _open = partial(gzip.open, mode='rt') if gzipped == True else open

        # iterate over the entries in input_colors
        with _open(input_colours[colour]) as handle:
            gff_record_list = []
            record_id = 1
            for record in SeqIO.parse(handle, "fasta"):
                # sort records based on first position
                GFF_entries[record_id].sort(key=lambda x: x[1][0][0])
                gff_record = record
                gff_record.features = []
                entry_ID = 1

                # remove duplicated entries
                seen = set()
                for entry in GFF_entries[record_id]:
                    start = entry[1][0][0] + 1
                    end = entry[1][0][1]
                    strand = 1 if entry[1][1] else -1

                    # ignore duplicates
                    ORF_position = (start, end, strand)
                    if ORF_position in seen:
                        continue
                    seen.add(ORF_position)

                    qualifiers = {
                        "source": "ggCaller:" + __version__,
                        "ID": isolate_names[colour] + "_" + str(entry[0]),
                        "inference": entry[2][0],
                        "score": entry[2][2],
                        "annotation": [entry[2][3].replace(",", " ").replace("|", "-")
                                           .replace("=", "-").replace("(", "").replace(")", "")
                                           .replace("[", "").replace("]", "").replace("'", "")
                                           .replace("<", "").replace(">", "").replace("~", "-")
                                           .replace(";", " ")],
                    }
                    feature = SeqFeature(
                        FeatureLocation(start, end, strand=strand),
                        type="CDS", qualifiers=qualifiers
                    )
                    gff_record.features.append(feature)
                    entry_ID += 1
                gff_record_list.append(gff_record)
                record_id += 1

        # write to GFF file
        gff_record_list = (x for x in gff_record_list)
        outfile = GFF_dir + isolate_names[colour] + ".gff"

        with open(outfile, "w") as out_handle:
            GFF.write(gff_record_list, out_handle, include_fasta=True)

    return