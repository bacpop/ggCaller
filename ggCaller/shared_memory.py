import sys
import numpy as np
import collections

try:
    from multiprocessing import Pool, shared_memory
    from multiprocessing.managers import SharedMemoryManager

    NumpyShared = collections.namedtuple('NumpyShared', ('name', 'shape', 'dtype'))
except ImportError as e:
    sys.stderr.write("This version of ggCaller requires python v3.8 or higher\n")
    sys.exit(1)


# generate shared memory array
def generate_shared_mem_array(in_array, smm):
    array_raw = smm.SharedMemory(size=in_array.nbytes)
    array_shared = np.ndarray(in_array.shape, dtype=in_array.dtype, buffer=array_raw.buf)
    array_shared[:] = in_array[:]
    array_shared_tup = NumpyShared(name=array_raw.name, shape=in_array.shape, dtype=in_array.dtype)
    return (array_shared, array_shared_tup)
