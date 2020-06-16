import sys
import subprocess
import run_profiler

def create_graph_file(nodes, path):
    graph_file = open(path, "w")
    graph_command = "python3 randomgraph.py -s {} --with-sequence".format(nodes)
    subprocess.call(graph_command, shell=True, stdout=graph_file)
    graph_file.close

def main(multipliers, result_file):
    file_handler = open(result_file, "w")
    nodes_previous_run = -1 # need to avoid computing duplication
    
    for multiplier in range(0, multipliers + 1):
        factor = 10 ** multiplier
        for number in [1, 2.5, 5, 7.5, 10]:            
            nodes = int(factor * number)
            if nodes == nodes_previous_run:
                continue
            path = "graph_nodes_{}.gfa".format(nodes)
            create_graph_file(nodes, path)
            result = run_profiler.run_profiler(path, "\n")
            file_handler.write(result)
            nodes_previous_run = nodes
    file_handler.close()

if __name__ == "__main__":
    multipliers = int(sys.argv[1])
    result_file = sys.argv[2]
    main(multipliers, result_file)
    subprocess.call("rm graph_nodes_*", shell=True)
