import sys
import os
import re
from operator import itemgetter

def read_gsm(n):
	f = open(n)
	by_query = {}
	q = None
	queries = []
	for l in f:
		l = l.rstrip("\n")
		if l.startswith("Query"):
			q = l.split(" ", 1)[1]
			queries.append(q)
			by_query.setdefault(q, [])
			continue
		ll = l.split()
		g = ll[0]
		pval1 = float(ll[1])
		pval2 = float(ll[2])
		if pval1<0.05 and pval2<0.05:
			by_query[q].append((g, pval1, pval2))

	f.close()
	return by_query, queries

def read_list(n):
	m = []
	f = open(n)
	for l in f:
		l = l.rstrip("\n")
		m.append(l)
	f.close()
	return m

if __name__=="__main__":
	EA = read_list("EA.list")
	AA = read_list("AA.list")
	#EA = read_list("ea.list.good")
	#AA = read_list("aa.list.good")
	EA_list = {}
	AA_list = {}
	queries = []
	for sample in EA:
		EA_list[sample], queries = read_gsm(sample)
	for sample in AA:
		AA_list[sample], queries = read_gsm(sample)

	#print(queries)
	#sys.exit(0)

	for ix,q in enumerate(queries):
		print("Query:", q)
		by_gene_EA = {}
		for sample in EA:
			for g, i, j in EA_list[sample][q]:
				by_gene_EA.setdefault(g, 0)
				by_gene_EA[g]+=1
		by_gene_AA = {}
		for sample in AA:
			for g, i, j in AA_list[sample][q]:
				by_gene_AA.setdefault(g, 0)
				by_gene_AA[g]+=1
		
		genes = set(list(by_gene_EA.keys()) + list(by_gene_AA.keys()))
		gs = []
		for g in genes:
			#print(q, g, by_gene_EA.get(g, 0), by_gene_AA.get(g, 0))
			count_EA = by_gene_EA.get(g, 0)
			count_AA = by_gene_AA.get(g, 0)
			size_ratio = len(AA) / len(EA)

			if count_EA==0:
				ratio = count_AA / len(AA) / (1/(len(EA)*2))
				score = ratio * (1.0 * size_ratio +count_AA) / 2
			else:
				ratio = count_AA/ len(AA)  / (count_EA / len(EA))
				score = ratio * (count_EA * size_ratio +count_AA) / 2

			#ratio>1.2
			if ratio>2.0 and count_AA/len(AA) >= 0.5:
			#if ratio<1.2 and ratio>1 and count_AA/len(AA) >= 0.5:
				gs.append((g, ratio, score, count_EA, count_AA))
		gs.sort(key=itemgetter(2), reverse=True)
		for g,c1,c2,c3,c4 in gs:
			print(g, c1, c2, c3, c4)
	
		#fw = open("q%d.shared2.more.txt" % (ix+1), "w")
		fw = open("M%d.AA.more.txt" % (ix+1), "w")
		for g,c1,c2,c3,c4 in gs:
			fw.write(g + "\n")
		fw.close()
