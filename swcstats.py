import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import networkx as nx

def swc_stats(filename, scale = 'mum', log=False):
    a = pd.read_csv(filename, sep=' ', header=None, comment='#')
    X = a.values
    if X.shape[1]>7:
        X = X[:, X.shape[1]-7:]
    G = nx.DiGraph()
    distance = 0
    surface_area = 0
    volume = 0
    for i in range(X.shape[0]):
        if X[i,6] != -1:
            G.add_node(i)
            parent = np.where(X[:,0] == X[i,6])[0][0]
            x_parent = X[parent,2:5]
            x = X[i,2:5]
            h = np.sqrt(np.sum(np.square(x_parent-x)))
            G.add_edge(parent,i,weight=h)
            distance += h
            r_parent = X[parent,5]
            r = X[i,5]
            surface_area += np.pi * (r + r_parent) * np.sqrt(np.square(r-r_parent)+np.square(h))
            volume += np.pi/3.*(r*r+r*r_parent+r_parent*r_parent)*h

    XX = X[:,2:5]
    w = np.abs(np.max(XX[:,0])-np.min(XX[:,0]))
    h = np.abs(np.max(XX[:,1])-np.min(XX[:,1]))
    d = np.abs(np.max(XX[:,2])-np.min(XX[:,2]))
    bifurcations = len(X[:,6])-len(np.unique(X[:,6]))
    max_euclidean_dist = np.max(pdist(XX))
    max_path_dist = nx.dag_longest_path_length(G)
    if log == True:
        print('Total Length: ', distance, scale)
        print('Total Surface Area: ', surface_area, scale+'^2')
        print('Total Volume: ', volume, scale+'^3')
        print('Maximum Euclidean Distance: ', max_euclidean_dist, scale)
        print('Width (Orientation Variant): ', w, scale)
        print('Height (Orientation Variant): ', h, scale)
        print('Depth (Orientation Variant): ', d, scale)
        print('Average Diameter: ', 2*np.mean(X[:,5]), scale)
        print('Number of Bifurcations:', bifurcations)
        print('Max Path Distance: ', max_path_dist, scale)
    results = {}
    results['Total Length'] = distance
    results['Total Surface Area'] = surface_area
    results['Total Volume'] = volume
    results['Maximum Euclidean Distance'] = max_euclidean_dist
    results['Width (Orientation Variant)'] = w
    results['Height (Orientation Variant)'] = h
    results['Depth (Orientation Variant)'] = d
    results['Average Diameter'] = 2*np.mean(X[:,5])
    results['Number of Bifurcations'] = bifurcations
    results['Max Path Distance'] = max_path_dist
    return results