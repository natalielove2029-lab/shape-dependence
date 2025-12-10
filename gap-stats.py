import numpy as np
import pandas as pd
import sys
import itertools
import glob
import os

L = sys.argv[1]
L = int(L)

kind = sys.argv[2]

###batch_max = sys.argv[3]
batch_max = 100 

###max_samples = sys.argv[4]
max_samples = sys.argv[3]
max_samples = int(max_samples)
batches = range(batch_max)

l = L//2

def cluster_gaps(cluster):
    n = np.zeros(l, dtype=np.float128)
    current_site = cluster[0]
    current_y = current_site // L
    first_site_in_row = current_site

    prev_site = current_site
    prev_y = current_y
    for ds in cluster[1:]:
        current_site = prev_site + ds
        current_y = current_site // L
        if current_y == prev_y:
            s = min(ds, L-ds)
            n[s - 1] += 1
        else:
            ds = prev_site - first_site_in_row
            first_site_in_row = current_site
            s = min(ds, L-ds)
            if s != 0:
                n[s - 1] += 1
        prev_site = current_site
        prev_y = current_y

    ds = prev_site - first_site_in_row
    s = min(ds, L-ds)
    if s != 0:
        n[s - 1] += 1

    return n

for b in batches:
    gaps = np.zeros(l, dtype=np.float128)
    samples = f'jonahSamples/ising/{kind}/{L}/{b}/*'
    count=0
    for file in glob.iglob(samples, recursive=True):
        if (count==max_samples):
            break
        with open(file, 'r') as f:
            fsize = os.path.getsize(file)
            if (fsize == 0):
                continue
            for raw_cluster in f:
                if (len(raw_cluster.strip()) > 0):
                    cluster = np.array(raw_cluster.split(), dtype=int)
                    gaps = gaps + cluster_gaps(cluster)
            count+=1
    output = f'jonahSamples/gap-stats/gaps_{L}_{kind}_{b}.txt'
    with open(output, 'w') as f:
        for num in gaps:
            f.write(f"{num}"+"\n")


