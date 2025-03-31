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
		avgexpr_h = float(ll[3])
		avgexpr_l = float(ll[4])
		a_mean = float(ll[5])
		a_std = float(ll[6])
		#if pval1<1e-10 and pval2<1e-10:
		if pval1<0.05 and pval2<0.05:
		#	by_query[q].append((g, pval1, pval2, avgexpr))
		#by_query[q].append((g, pval1, pval2, avgexpr_h/avgexpr_l))
			by_query[q].append((g, pval1, pval2, (avgexpr_h - a_mean)/a_std))
			#by_query[q].append((g, pval1, pval2, avgexpr_h - a_mean))
			#by_query[q].append((g, pval1, pval2, avgexpr_h / a_std))

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
	EA_list = {}
	AA_list = {}
	queries = []
	for sample in EA:
		EA_list[sample], queries = read_gsm(sample)
	for sample in AA:
		AA_list[sample], queries = read_gsm(sample)

	#print(queries)
	#sys.exit(0)

	#genes_x = read_list(sys.argv[1])
	#genes_x = read_list("selected.by.expr.q%s.EA.more.txt" % sys.argv[1])
	genes_x = read_list("selected.by.expr.q%s.EA.more.with.neutro.txt" % sys.argv[1])

	qid = int(sys.argv[1])
	by_gene = {}
	for ix,q in enumerate(queries):
		#print("Query:", q)

		by_gene_EA = {}
		for sample in EA:
			for g, i, j, a in EA_list[sample][q]:
				by_gene_EA.setdefault(g, {})
				by_gene_EA[g].setdefault(sample, 0)
				by_gene_EA[g][sample] = a
		by_gene_AA = {}
		for sample in AA:
			for g, i, j, a in AA_list[sample][q]:
				by_gene_AA.setdefault(g, {})
				by_gene_AA[g].setdefault(sample, 0)
				by_gene_AA[g][sample] = a

		if ix == qid-1:
			for g in genes_x:
				sam_EA = []
				sam_AA = []
				for s in EA:
					if g not in by_gene_EA or s not in by_gene_EA[g]:
						sam_EA.append(0)
					else:
						sam_EA.append(by_gene_EA[g][s])
				for s in AA:
					if g not in by_gene_AA or s not in by_gene_AA[g]:
						sam_AA.append(0)
					else:
						sam_AA.append(by_gene_AA[g][s])

				#sam_EA = [by_gene_EA[g][s] for s in EA]
				#sam_AA = [by_gene_AA[g][s] for s in AA]
				sys.stdout.write(g + "\t" + "\t".join(["%.2f" % x for x in sam_AA]) + "\t" + "\t".join(["%.2f" % x for x in sam_EA]) + "\n")
		
	sys.exit(0)
		
	'''
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

			if ratio>1.2 and count_AA/len(AA) >= 0.5:
			#if ratio<1.2 and ratio>1 and count_AA/len(AA) >= 0.5:
				gs.append((g, ratio, score, count_EA, count_AA))
		gs.sort(key=itemgetter(2), reverse=True)
		for g,c1,c2,c3,c4 in gs:
			if g in genes_x:
				#print(g, c1, c2, c3, c4)
				by_gene.setdefault(g, {})
				by_gene[g][ix] = c1
	'''
	
	'''
		fw = open("q%d.shared2.more.txt" % (ix+1), "w")
		for g,c1,c2,c3,c4 in gs:
			fw.write(g + "\n")
		fw.close()
	'''

	for g in genes_x:
		if g not in by_gene: continue
		sys.stdout.write(g + "\t")
		for ix in range(9):	
			sys.stdout.write(str(by_gene[g].setdefault(ix, 0)) + "\t")
		sys.stdout.write("\n")
