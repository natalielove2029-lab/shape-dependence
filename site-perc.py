import sys
import numpy as np
import networkx as nx

size = sys.argv[1]
size = int(size)
num = sys.argv[2]
p = 0.592746
G = nx.grid_2d_graph(size, size, periodic=True)
G = nx.convert_node_labels_to_integers(G, ordering='sorted')
for u in range(size**2):
    if (np.random.rand() > p):
        G.remove_node(u)

outfile = f'./perc-site-{size}-{num}.txt'

for i in range(size**2):
    if G.has_node(i):
        cluster = nx.node_connected_component(G, i)
        G.remove_nodes_from(cluster)
        sites = sorted(cluster)
        cluster_id = sites[0]
        prev = cluster_id
        if (len(sites) > 1):
            with open(outfile, 'a+') as f:
                f.write(f'{cluster_id}' + '\t')
                for site in sites[1::]:
                    ds = site - prev
                    f.write(f'{ds}' + '\t')
                    prev = site
                f.write('\n')