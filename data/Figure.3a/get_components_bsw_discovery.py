import sys
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
def read(n):
	f = open(n)
	f.readline()
	edges_AA = []
	edges_EA = []
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		aa_score = float(ll[1])
		ea_score = float(ll[2])
		edge = ll[0].split("--")
		edge[0] = int(edge[0].split()[0])
		edge[1] = int(edge[1].split()[0])
		if aa_score>ea_score:
			edges_AA.append((edge[0], edge[1]))
		else:
			edges_EA.append((edge[0], edge[1]))
	f.close()
	return edges_AA, edges_EA

if __name__=="__main__":
	e_AA, e_EA = read("bsw.discovery")
	
	G = nx.Graph()
	G.add_edges_from(e_AA)
	G.add_edges_from([(20,19),]) #joined, due to influence from other datasets (BSW-validation, and Roswell)

	components = list(nx.connected_components(G))
	print("AA modules")
	for i, component in enumerate(components):
		print(f"Component {i+1}: {component}")

	node_color = {}
	for c in [4, 12, 14, 20, 1, 19, 13, 9]:
		node_color[c] = "red"
	for c in [3, 5, 6, 8, 16, 18]:
		node_color[c] = "blue"
	node_colors = [node_color[c] for c in G.nodes()]

	# Create side-by-side plots
	fig, axes = plt.subplots(2, 2, figsize=(14, 12))
	#plt.figure(figsize=(8, 6))
	axes = axes.flatten()

	#plt.figure(figsize=(8, 6))
	pos = nx.spring_layout(G, seed=42)  # layout for better visualization
	dash_edges = [(20,19),]
	normal_edges = [e for e in G.edges if e not in dash_edges and (e[1], e[0]) not in dash_edges]

	nx.draw_networkx_nodes(G, pos, ax=axes[0], node_color=node_colors, node_size=400)
	nx.draw_networkx_labels(G, pos, ax=axes[0])

	nx.draw_networkx_edges(G, pos, ax=axes[0], edgelist=normal_edges, edge_color='gray')
	nx.draw_networkx_edges(G, pos, ax=axes[0], edgelist=dash_edges, edge_color='red', style='dotted', width=2)

	legend_handles = []
	patch = Patch(color="red", label=f"BA Community 1")
	legend_handles.append(patch)
	patch = Patch(color="blue", label=f"BA Community 2")
	legend_handles.append(patch)
	axes[0].legend(handles=legend_handles, title="Legend:", loc='upper left')
	axes[0].set_title("BA specific communities")

	# Get adjacency matrix (as a NumPy array)
	nodes = [4, 12, 14, 20, 1, 19, 13]
	Gsub = G.subgraph(nodes).copy()
	Gsub.remove_edges_from([(20,19),]) #remove artificial edges
	adj_matrix = nx.to_numpy_array(Gsub, nodelist=nodes)
	df_adj = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)
	sns.heatmap(df_adj, cmap="Greens", cbar=False, linewidths=0.5, linecolor='gray', square=True, ax=axes[2])
	axes[2].set_title("BA-community 1 Adjacency Matrix Heatmap")



	#=========================================================================
	# Add legend to second plot (or wherever you prefer)
	axes[0].legend(handles=legend_handles, title="Legend:", loc='upper left')
	axes[0].set_title("BA specific communities")

	G = nx.Graph()
	G.add_edges_from(e_EA)
	G.add_edges_from([(2,8),]) #joined, due to influence from other datasets (BSW-validation, and Roswell)
	components = list(nx.connected_components(G))
	print("EA modules")
	for i, component in enumerate(components):
		print(f"Component {i+1}: {component}")

	node_color = {}
	for c in [11, 15, 9, 7, 12, 10, 14, 6]:
		node_color[c] = "red"
	for c in [16, 17, 2, 8, 19]:
		node_color[c] = "blue"
	for c in [1, 13]:
		node_color[c] = "lightblue"	
	node_colors = [node_color[c] for c in G.nodes()]

	pos = nx.spring_layout(G, seed=4)  # layout for better visualization
	dash_edges = [(2,8),]
	normal_edges = [e for e in G.edges if e not in dash_edges and (e[1], e[0]) not in dash_edges]
	nx.draw_networkx_nodes(G, pos, ax=axes[1], node_color=node_colors, node_size=400)
	nx.draw_networkx_labels(G, pos, ax=axes[1])
	nx.draw_networkx_edges(G, pos, ax=axes[1], edgelist=normal_edges, edge_color='gray')
	nx.draw_networkx_edges(G, pos, ax=axes[1], edgelist=dash_edges, edge_color='red', style='dashed', width=2)

	legend_handles = []
	patch = Patch(color="red", label=f"WA Community 1")
	legend_handles.append(patch)
	patch = Patch(color="blue", label=f"WA Community 2")
	legend_handles.append(patch)
	axes[1].legend(handles=legend_handles, title="Legend:", loc='upper right')
	axes[1].set_title("WA specific communities")

	# Get adjacency matrix (as a NumPy array)
	nodes = [11, 15, 9, 7, 12, 10, 14, 6]
	Gsub = G.subgraph(nodes).copy()
	adj_matrix = nx.to_numpy_array(Gsub, nodelist=nodes)
	df_adj = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)
	sns.heatmap(df_adj, cmap="Greens", cbar=False, linewidths=0.5, linecolor='gray', square=True, ax=axes[3])
	axes[3].set_title("WA-community 1 Adjacency Matrix Heatmap")
	
	#nx.draw(G, pos, with_labels=True, node_color="red", cmap=plt.cm.Set3, node_size=400, edge_color='gray', ax=axes[1])
	plt.tight_layout()
	plt.show()
