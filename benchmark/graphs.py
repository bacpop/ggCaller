"""
Python script to create graphs of the benchmark
"""

from collections import OrderedDict
import matplotlib.pyplot as plt

BENCHMARK_NODE = OrderedDict()
with open('results/cmp_nodes.txt') as node_method_file:
    ROWS = node_method_file.readlines()
    i = 0
    for ROW in ROWS:
        i += 1
        ROW = ROW.split('\t')
        if not BENCHMARK_NODE.get(int(ROW[6])):
            BENCHMARK_NODE[int(ROW[6])] = {i : {'overlap_consistency' : float(ROW[2])}}
        else:
            BENCHMARK_NODE[int(ROW[6])][i] = {'overlap_consistency' : float(ROW[2])}

        BENCHMARK_NODE[int(ROW[6])][i]['node_compression'] = float(ROW[3])
    node_method_file.close()

ELEMENT_NODE = []
CHECK_NODE =[]
CMP_NODE = []
for e in sorted(BENCHMARK_NODE.keys()):
    for value in BENCHMARK_NODE.get(e):
        ELEMENT_NODE.append(e)
        CHECK_NODE.append(BENCHMARK_NODE.get(e).get(value).get('overlap_consistency'))
        CMP_NODE.append(BENCHMARK_NODE.get(e).get(value).get('node_compression'))


#plt.plot(ELEMENT_NODE, CHECK_NODE, label='Overlap consistency')

plt.plot(ELEMENT_NODE, CMP_NODE, label='Compression by node')

import numpy as np

x = np.linspace(0, 2, 100)
sec = np.power(x, 2)
terz = np.power(x, 3)
e = np.exp(x)
e2 = np.exp2(x)

#plt.plot(x,x, label='Lin')
#plt.plot(x,sec, label='Sec')
#plt.plot(x,terz, label='Terz')
#plt.plot(x,e, label='Exp')
#plt.plot(x,e2, label='Exp2')


plt.xlabel('number nodes + edges')
plt.ylabel('seconds')
plt.legend()
plt.savefig("node_compression.png", dpi=500)
plt.show()
