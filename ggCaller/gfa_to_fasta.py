def gfa_to_fasta(in_file):
    """Converts .gfa file to .fasta, using same file name.

    Args:
        in_file (str)
            path for .gfa file to convert.
    """
    base = os.path.splitext(in_file)[0]
    out_file = base + ".fasta"
    with open(in_file, "r") as f, open(out_file, "w") as o:
        for line in f:
            parsed_line = re.split(r'\t+', line.rstrip('\t'))
            if parsed_line[0] == "S":
                    o.write(">" + str(parsed_line[1]) + "\n" + str(parsed_line[2]) + "\n")

if __name__ == '__main__':
    import sys
    import os
    import re

    graph_file = sys.argv[1]
    gfa_to_fasta(graph_file)