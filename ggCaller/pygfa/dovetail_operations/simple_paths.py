from pygfa.algorithms.simple_paths import all_simple_paths

def dovetails_all_simple_paths(gfa_, source, target, edges=False, keys=False, cutoff=None):
    return all_simple_paths(\
                            gfa_, \
                            source, \
                            target, \
                            gfa_.dovetails_iter, \
                            edges=edges, \
                            keys=keys, \
                            cutoff=cutoff) # argument for gfa_dovetails_iter
