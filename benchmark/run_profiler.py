import time
import sys
from pympler import asizeof
sys.path.insert(1, '../')
import pygfa


def timeit(function):
    def timed(*args, **kwargs):
        kw = kwargs.copy()
        if "log_data" in kwargs:
            kw.pop("log_data")
        ts = time.time()
        result = function(*args, **kw)
        te = time.time()
        time_ = "{0:f}".format(te-ts)
        if "log_data" in kwargs:
            kwargs["log_data"].append(time_)
        else:
            print(time_)
        return result
    return timed

@timeit
def load_graph(file_path):
    gfa_ = pygfa.gfa.GFA.from_file(file_path)
    return gfa_

@timeit
def compute_elements(gfa_):
    nodes = gfa_.nodes()
    edges = gfa_.edges()
    x = [x for x in range(1, 2**10)]
    return len(nodes), len(edges)

@timeit
def compute_connected_components(gfa_):
    conn_components = list(pygfa.nodes_connected_components(gfa_))
    dov_conn_components = list(pygfa.dovetails_nodes_connected_components(gfa_))
    return len(conn_components), len(dov_conn_components)

@timeit
def compute_linear_paths(gfa_):
    linear_paths = list(pygfa.dovetails_linear_paths(gfa_))
    return len(linear_paths)

@timeit
def compute_overlap_consistency(gfa_):
    edges_no_consistency, edges_no_calculate = pygfa.gfa.GFA.overlap_consistency(gfa_)
    return len(edges_no_consistency), len(edges_no_calculate)

@timeit
def compute_compression_by_nodes(gfa_):
    before_n_edges = len(gfa_.edges())
    pygfa.gfa.GFA.compression(gfa_)
    after_n_edges = len(gfa_.edges())
    return (before_n_edges, after_n_edges)

@timeit
def compute_compression_by_edges(gfa_):
    before_n_edges = len(gfa_.edges())
    pygfa.gfa.GFA.compression(gfa_, 'by_edges')
    after_n_edges = len(gfa_.edges())
    return (before_n_edges, after_n_edges)

def run_profiler(file_, end=""):
    data = []
    gfa_ = load_graph(file_, log_data=data)
    nodes, edges = compute_elements(gfa_, log_data=data)
    cc, dov_cc = compute_connected_components(gfa_, log_data=data)
    lin_paths = compute_linear_paths(gfa_, log_data=data)
    data.extend([nodes, edges, cc, dov_cc, lin_paths, asizeof.asizeof(gfa_)])
    return str.join("\t", [str(x) for x in data]) + end

def run_profiler_graph_operation(file_, end="", type_compression='nodes'):
    data = []
    gfa_ = load_graph(file_, log_data=data)
    nodes, edges = compute_elements(gfa_, log_data=data)
    elements = nodes + edges
    ovr_cons = compute_overlap_consistency(gfa_, log_data=data)
    compr = None
    if type_compression == 'edges':
        compr = compute_compression_by_edges(gfa_, log_data=data)
    else:
        compr = compute_compression_by_nodes(gfa_, log_data=data)
    data.extend([nodes, edges, elements, ovr_cons, compr])
    return str.join("\t", [str(x) for x in data]) + end

if __name__ == "__main__":
    print(run_profiler(sys.argv[1]))
