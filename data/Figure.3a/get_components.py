import sys
import networkx as nx
import matplotlib.pyplot as plt
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
	fig, axes = plt.subplots(1, 2, figsize=(14, 6))
	plt.figure(figsize=(8, 6))
	pos = nx.spring_layout(G, seed=42)  # layout for better visualization
	nx.draw(G, pos, with_labels=True, node_color="red", cmap=plt.cm.Set3, node_size=400, edge_color='gray', ax=axes[0])
	#plt.title("Connected Components in Graph")
	#plt.show()

	G = nx.Graph()
	G.add_edges_from(e_EA)
	components = list(nx.connected_components(G))
	print("EA modules")
	for i, component in enumerate(components):
		print(f"Component {i+1}: {component}")
	
	pos = nx.spring_layout(G, seed=4)  # layout for better visualization
	nx.draw(G, pos, with_labels=True, node_color="red", cmap=plt.cm.Set3, node_size=400, edge_color='gray', ax=axes[1])
	plt.tight_layout()
	plt.show()
