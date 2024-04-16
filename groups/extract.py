import sys
import os
import re

f = open("full.matrix.txt")
m = {}
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	t1 = ll[0]
	m.setdefault(t1, [])
	m[t1].append(l)
f.close()
for k in m:
	k2 = k.replace("/", " ")
	fw = open("%s.txt" % k2, "w")
	for l in m[k]:
		fw.write(l + "\n")
	fw.close()
