import os
import numpy as np
import pandas as pd
import networkx as nx

## uncomment for plotting (e.g., investigating networks)
# import matplotlib
# matplotlib.use("TkAgg")
# from matplotlib import pyplot as plt

# def drawit(g, **kwargs):
#     nx.draw(g)
#     plt.show()

# local directory for the ASNR repository
# https://github.com/bansallab/asnr
topdir = "../asnr/Networks/Mammalia/"
netdirs = os.listdir(topdir)

results = []

# "baboon_franz_gooming_group_12" produces an NA value below. Here is the network:
# EdgeDataView([('4', '7', {'weight': 1}), ('5', '6', {'weight': 1})])
# This network is disconnected, so NA behavior is appropriate

# investigate every network in topdir, calculating (b/c)* as able
for subdir in netdirs:
    print(subdir)

    # Skip bipartite networks
    if "bipartite" in subdir:
        continue

    # collect the appropriate networks
    train_path = os.path.join(topdir, subdir)
    graphml_list = [
        os.path.join(train_path, i) for i in os.listdir(train_path) if i.endswith(".graphml")
    ]

    # and check each one
    for iteration in range(len(graphml_list)):
        # Load the graph
        G = nx.read_graphml(graphml_list[iteration])
        # take a largest connected if it exists
        try:
            largest = max(nx.connected_components(G), key = len)
            G = G.subgraph(list(largest))
        except:
            continue
        # skipping networks that are too large for the below calculations
        if G.number_of_nodes() > 100:
            continue

        # the next several lines (to `Threshold = t2/(t3 - t1)`) implement the
        # procedure of Allen et al 2017, "Evolutionary dynamics on any population structure"
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

        # Calculate and store the threshold and relevant network metrics
        result = {
            "threshold": Threshold,
            "network": graphml_list[iteration],
            "n_nodes": G.number_of_nodes(),
            "mean_degree": sum([d[1] for d in G.degree])/G.number_of_nodes(),
            "avg_clustering": nx.average_clustering(G),
            "mean_strength": sum(d[1] for d in G.degree(weight = "weight"))/G.number_of_nodes() if nx.is_weighted(G) else np.nan,
            "avg_weighted_clustering": nx.average_clustering(G, weight = "weight") if nx.is_weighted(G) else np.nan
        }
        results.append(result)

# make the above into a data frame for export
cleaned = []
for i in range(len(results)):
    bc = results[i]["threshold"][0, 0]
    N = results[i]["n_nodes"]
    kmean = results[i]["mean_degree"]
    clust = results[i]["avg_clustering"]
    smean = results[i]["mean_strength"]
    wclust = results[i]["avg_weighted_clustering"]
    net = results[i]["network"]
    cleaned.append({"bc": bc, "N": N, "kmean": kmean, "clust": clust,
                    "smean": smean, "wclust": wclust,
                    "net": net})

df = pd.DataFrame(cleaned)
df.to_csv("network-data.csv")
print("Done.")
