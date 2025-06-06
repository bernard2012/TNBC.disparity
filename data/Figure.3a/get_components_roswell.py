import sys
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.patches import Wedge
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
	e_AA, e_EA = read(sys.argv[1])
	G = nx.Graph()
	G.add_edges_from(e_AA)
	components = list(nx.connected_components(G))
	print("AA modules")
	for i, component in enumerate(components):
		print(f"Component {i+1}: {component}")

	# Create side-by-side plots
	fig, axes = plt.subplots(2, 2, figsize=(14, 12))
	axes = axes.flatten()
	pos = nx.spring_layout(G, seed=42)  # layout for better visualization

	node_color = {}
	for c in [2, 6, 9, 14, 15, 8, 19]:
		node_color[c] = "red"
	for c in [4, 20, 16, 11, 10]:
		node_color[c] = "blue"
	for c in [5, 7, 17]:
		node_color[c] = "green"
	node_colors = [node_color[c] for c in G.nodes()]

	nx.draw(G, pos, with_labels=True, node_color=node_colors, cmap=plt.cm.Set3, node_size=400, edge_color='gray', ax=axes[0])
	#plt.title("Connected Components in Graph")
	#plt.show()
	# Build legend handles
	legend_handles = []
	patch = Patch(color="red", label=f"BA Community 1")
	legend_handles.append(patch)
	patch = Patch(color="blue", label=f"BA Community 2")
	legend_handles.append(patch)
	patch = Patch(color="green", label=f"BA Community 3")
	legend_handles.append(patch)

	# Add legend to second plot (or wherever you prefer)
	axes[0].legend(handles=legend_handles, title="Legend:", loc='upper left')
	axes[0].set_title("BA specific communities")

	# Get adjacency matrix (as a NumPy array)
	nodes = [15, 9, 6, 8, 2, 14, 19]
	Gsub = G.subgraph(nodes).copy()
	#Gsub.remove_edges_from([(20,19),]) #remove artificial edges
	adj_matrix = nx.to_numpy_array(Gsub, nodelist=nodes)
	df_adj = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)
	sns.heatmap(df_adj, cmap="Greens", cbar=False, linewidths=0.5, linecolor='gray', square=True, ax=axes[2])
	axes[2].set_title("BA-community 1 Adjacency Matrix Heatmap")

	#===================================================
	G = nx.Graph()
	G.add_edges_from(e_EA)
	components = list(nx.connected_components(G))
	print("EA modules")
	for i, component in enumerate(components):
		print(f"Component {i+1}: {component}")

	special_node = 9

	node_color = {}
	for c in [14, 15, 16, 18, 10, 17, 1, 9]:
		node_color[c] = "red"
	for c in [11, 13, 19]:
		node_color[c] = "green"
	for c in [2, 5, 6, 7, 8]:
		node_color[c] = "blue"
	for c in [3]:
		node_color[c] = "lightblue"

	pos = nx.spring_layout(G, seed=4)  # layout for better visualization
	default_nodes = [n for n in G.nodes if n != special_node]
	node_colors = [node_color[c] for c in default_nodes]
	nx.draw_networkx_edges(G, pos, edge_color='gray', ax=axes[1])
	nx.draw_networkx_nodes(G, pos, nodelist=default_nodes, node_color=node_colors, node_size=400, ax=axes[1])
	nx.draw_networkx_labels(G, pos, ax=axes[1])

	radius = 0.05
	x,y = pos[special_node]
	wedge1 = Wedge(center=(x, y), r=radius, theta1=90, theta2=270, facecolor='red', edgecolor='black')
	wedge2 = Wedge(center=(x, y), r=radius, theta1=270, theta2=90, facecolor='green', edgecolor='black')

	axes[1].add_patch(wedge1)
	axes[1].add_patch(wedge2)

	axes[1].text(x, y, str(special_node), fontsize=12, ha='center', va='center', color='white')
	axes[1].set_title("WA specific communities")

	legend_handles = []
	patch = Patch(color="red", label=f"WA Community 1")
	legend_handles.append(patch)
	patch = Patch(color="blue", label=f"WA Community 2")
	legend_handles.append(patch)
	patch = Patch(color="green", label=f"WA Community 3")
	legend_handles.append(patch)

	# Add legend to second plot (or wherever you prefer)
	axes[1].legend(handles=legend_handles, title="Legend:", loc='lower right')


	# Get adjacency matrix (as a NumPy array)
	nodes = [1, 17, 10, 9, 16, 18, 15, 14]
	Gsub = G.subgraph(nodes).copy()
	adj_matrix = nx.to_numpy_array(Gsub, nodelist=nodes)
	df_adj = pd.DataFrame(adj_matrix, index=nodes, columns=nodes)
	sns.heatmap(df_adj, cmap="Greens", cbar=False, linewidths=0.5, linecolor='gray', square=True, ax=axes[3])
	axes[3].set_title("WA-community 1 Adjacency Matrix Heatmap")


	#nx.draw(G, pos, with_labels=True, node_color="red", cmap=plt.cm.Set3, node_size=400, edge_color='gray', ax=axes[1])
	plt.tight_layout()
	plt.show()
