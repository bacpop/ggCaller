import sys
import subprocess
import run_profiler

def main(result_file, in_file, type_compression='nodes'):
    
    with open(in_file) as graph_file:
        list_graph = graph_file.readlines()
        graph_file.close()

    file_handler = open(result_file, "w")

    for file in list_graph:
        file = file.rstrip('\n')
        print(file)
        result = run_profiler.run_profiler_graph_operation(file, "\n", type_compression)
        file_handler.write(result)

    file_handler.close()

if __name__ == "__main__":
    result_file = sys.argv[1]
    graph_file = sys.argv[2]
    type_compression = sys.argv[3]
    main(result_file, graph_file, type_compression)
