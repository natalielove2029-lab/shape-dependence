import numpy as np
import pandas as pd
import sys
import os

### FUNCTIONS ###
#Note: 'square' is equivalent to sheared square, linear size l = L/2. All subsystems are represented as a 1d array of their cluster counts.
"""
Direct counting:
Given the location of the upper-left corner of a square
Outputs an array of number of contained sites of each cluster
"""
def square_inside_sites(pos0):
    x0 = pos0 % L
    y0 = pos0 // L
    
    pos = pos0
    cluster_counts = np.zeros(cn, dtype=int)
    
    for j in range(l):
        y = (y0 + j) % L
        for i in range(l):
            x = (x0 + i) % L
            pos = x + (y*L)
            cluster = sites[pos]
            if (cluster >= 0):
                cluster_counts[cluster] = cluster_counts[cluster] + 1
        x0 = (x0 + slant) % L
    return cluster_counts

"""
Fast slab:
Given 2 non-overlapping squares (square1, square2)
Outputs their union
"""
def slab_inside_sites(square1, square2):
    slab = np.zeros(cn, dtype=int)
    for i in range(cn):
        slab[i] = square1[i] + square2[i]
    return slab

"""
Update counting:
Given the location of the upper-left corner of a square, and the square located 1 to the left
Outputs an array of number of contained sites of each cluster
"""
def square_update(pos0, prev):
    cluster_counts = np.zeros(cn, dtype=int)
    for i in range(cn):
        cluster_counts[i] = prev[i]
    
    x0 = pos0 % L
    y0 = pos0 // L

    for j in range(l):
        x1 = (x0 - 1) % L
        y = (y0 + j) % L
        pos1 = x1 + y * L
        cluster = sites[pos1]
        if cluster >= 0:
            cluster_counts[cluster] = cluster_counts[cluster] - 1
            
        x2 = (x1 + l) % L
        pos2 = x2 + y * L
        cluster = sites[pos2]
        if cluster >= 0:
            cluster_counts[cluster] = cluster_counts[cluster] + 1

        x0 = (x0 + slant) % L
    return cluster_counts    

"""
Entanglement entropy:
Given a slab or square (subsystem)
Outputs entanglement entropy of that subsystem
"""
def entanglement_entropy(subsystem):
    s = 0
    for n in range(cn):
        if (0 < subsystem[n] < cluster_sizes[n]):
            s = s + 1
    return s

"""
One config. corner contribution:
Given 4 squares (sq1, sq2, sq3, sq4) which partition the sample according to the geom. method
Outputs the corner contribution of those 4 squares
"""
def corner_contrib(sq1, sq2, sq3, sq4):
    sl1 = slab_inside_sites(sq1, sq2)
    sl2 = slab_inside_sites(sq1, sq3)
    s_sq = entanglement_entropy(sq1) + entanglement_entropy(sq2) + entanglement_entropy(sq3) + entanglement_entropy(sq4)
    s_sl = 2 * (entanglement_entropy(sl1) + entanglement_entropy(sl2))
    s_corner = (s_sl - s_sq) / 4
    return s_corner

"""
Main program:
Outputs position-averaged S_cr
"""
def main():
    p1 = 0
    p2 = l
    x3 = (l * slant) % L
    p3 = x3 + (l * L)
    x4 = (x3 + l) % L
    p4 = x4 + (l * L)
    
    s1 = square_inside_sites(p1)
    s2 = square_inside_sites(p2)
    s3 = square_inside_sites(p3)
    s4 = square_inside_sites(p4)
    
    s_cr = np.zeros(l**2)
    
    for j in range(l):
        p1old = p1
        p2old = p2
        p3old = p3
        p4old = p4
        for i in range(l):
            s_cr[i + j * l] = corner_contrib(s1, s2, s3, s4)
            if (i == l - 1):
                p1 = p1old + L
                p2 = p2old + L
                p3 = p3old + L
                p4 = p4old + L
                s1 = square_inside_sites(p1)
                s2 = square_inside_sites(p2)
                s3 = square_inside_sites(p3)
                s4 = square_inside_sites(p4)
            else:
                p1 = p1 + 1
                p2 = p2 + 1
                p3 = p3 + 1
                p4 = p4 + 1
                s1 = square_update(p1, s1)
                s2 = square_update(p2, s2)
                s3 = square_update(p3, s3)
                s4 = square_update(p4, s4)
        
    s_cr_avg = sum(s_cr)/len(s_cr)
    return s_cr_avg

### SYSTEM ARGUMENTS ###
"""
file: filepath of sample
L: linear size of sample
disorder: disorder type (for internal file management)
slant: determines the shear angle, so that tan(γ) = 1/slant
batch: batch no of samples (for internal file management)
"""
file = sys.argv[1]
L = sys.argv[2]
disorder = sys.argv[3]
max_slant = sys.argv[4]

L = int(L)
if (L % 2 == 1):
    raise Exception("Invalid sample size")
max_slant = int(max_slant)
l = L//2 

### PREPROCESSING ###
"""
Data preprocessing: check if the sample is empty (then corner contribution is 0)
Otherwise, read the data in line-by-line; each line of raw data is a cluster.
Make dictionary keying cluster id to array index.
Make list of total cluster sizes.
Sample represented as L x L array of which cluster (by index) each site belongs to. -1 represents no cluster.
All explicit positions are in the form (x + y × L), where (x, y) are counted from the upper-left starting at 0.
Periodic BC
"""
fsize = os.path.getsize(file)

clust_ids = {}
ind = 0
cluster_sizes = np.zeros(0)
sites = np.full(L**2, -1, dtype=int)
cn = 0

with open(file, 'r+') as f:
    for raw_cluster in f:
        cluster = np.array(raw_cluster.split(), dtype=int)
        id = cluster[0]
        clust_ids[id] = cn
        size = len(cluster)
        cluster_sizes = np.append(cluster_sizes, size)
        site = id
        sites[site] = cn
        for i in range(1, len(cluster)):
            site = site + cluster[i]
            sites[site] = cn
        cn = cn + 1

### FUNCTION CALL ###

for slant in range(max_slant):
    if (fsize==0):
        scr = 0.0
    else:
        scr = main()
    outfile = f'scr-{L}-{disorder}-{slant}.txt'
    with open(outfile, 'a') as f:
        f.write(f'{scr}' + '\n')
