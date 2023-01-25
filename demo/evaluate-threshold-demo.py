import numpy as np
import networkx as nx

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

def eval_threshold(G):
    '''
    Returns a scalar (class 'float') value for the cooperation threshold (b/c)* of the (possibly weighted, possibly directed) network `G`.
    '''
    A = nx.to_numpy_matrix(G)
    D = np.sum(A, axis = 1)
    D = np.diagflat(D)
        
    P = np.dot(np.linalg.inv(D), A)
    P2 = np.dot(P, P)
    P3 = np.dot(P2, P)
        
    Pi = np.zeros(G.number_of_nodes())
    W = A.sum()
    for i in range(G.number_of_nodes()):
        Pi[i] = D[i, i]/W
        
    N = G.number_of_nodes()
    Mat = np.zeros((N**2, N**2))
    Mat = np.asmatrix(Mat)
    for i in range(N):
        for j in range(N):
            if i != j:
                Mat[N*i+j:N*i+j+1, N*j:N*(j+1)] = P[i:i+1]
                Mat[N*i+j:N*i+j+1, N*i:N*(i+1)] = P[j:j+1]
                Mat[N*i+j, N*j+i] = 0
                Mat[N*i+j, N*j+j] = 0
                Mat[N*i+j, N*i+i] = 0
                Mat[N*i+j, N*i+j] = 0
    Coef = np.asmatrix(np.identity(N**2)) - Mat/2
                    
    b = np.full(G.number_of_nodes()**2, 1)
        
    for i in range(G.number_of_nodes()**2):
        if i//G.number_of_nodes() == i%G.number_of_nodes():
            b[i] = 0

    x = np.linalg.solve(Coef, b)

    Tau = np.zeros((N, N))
    Tau = np.asmatrix(Tau)
    for i in range(N**2):
        Tau[i//N, i%N] = x[i]

    t1 = np.dot(Pi, np.sum(np.multiply(P, Tau), axis = 1))
    t2 = np.dot(Pi, np.sum(np.multiply(P2, Tau), axis = 1))
    t3 = np.dot(Pi, np.sum(np.multiply(P3, Tau), axis = 1))
    Threshold = t2/(t3 - t1)

    return float(Threshold)
    
def get_gcc(G):
    '''
    Return the largest connected component of `G`.
    '''
    largest = max(nx.connected_components(G), key = len)
    
    return G.subgraph(list(largest))

# The grooming network of a group of Cercopithecus campbelli, a typical catarrhine primate
# https://en.wikipedia.org/wiki/Campbell's_mona_monkey
# https://animaldiversity.org/accounts/Cercopithecus_campbelli/
# https://github.com/bansallab/asnr/tree/master/Networks/Mammalia/cercopithecus_campbelli_grooming_griffin
graphfile = "weighted_primate_matrix_5.graphml"
g = get_gcc(nx.read_graphml(graphfile))
k = [x[1] for x in g.degree]
w = [x[2]["weight"] for x in g.edges(data = True)]
pos = nx.spring_layout(g)#kamada_kawai_layout(g)
netstats = {
    "N": g.number_of_nodes(),
    "mean_degree": round(sum([d[1] for d in g.degree])/g.number_of_nodes(), 3),
    "avg_clustering": round(nx.average_clustering(g), 3),
    "mean_strength": round(sum(d[1] for d in g.degree(weight = "weight"))/g.number_of_nodes(), 3),
    "avg_weighted_clustering": round(nx.average_clustering(g, weight = "weight"), 3),
    "(b/c)*": round(eval_threshold(g), 3)
}
print(netstats)

fig, ax = plt.subplots()
nx.draw(
    g, with_labels = False, pos = pos, ax = ax,
    node_size = [50*x for x in k], node_color = "goldenrod",
    width = w, edge_color = "#36454f"
)
fig.show()
