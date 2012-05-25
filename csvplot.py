import dumpcnt
from pprint import pprint 
import re
from collections import OrderedDict
import csv
import sys

def transpose(d):
    td = OrderedDict([])
    for k1, v in d.items():
        for k2, v2 in v.items():
            if not k2 in td:
                td[k2] = OrderedDict([])
            if not k1 in td[k2]:
                td[k2][k1] = v2

    return td
                



f = csv.reader(open(sys.argv[1]))
mnlist = f.next()[1:]
data = []
for l in f:
    data.append(l)
d = OrderedDict([])
for mn in data:
    if not mn[0] in d:
        d[mn[0]] = []
for l in data:
    if l:
        d[l[0]].append([float(x) for x in l[1:]])
d1=OrderedDict()

for k in d:
    l = len(d[k][0])
    tmp = [[] for i in range(l)]
    for i in range(l):
        for x in d[k]:
            tmp[i].append(x[i])
    d1[k] = tmp
d = d1
pprint(d)
gd = OrderedDict()
for k in d:
    gd[k] = OrderedDict(zip(mnlist, [{"mean":dumpcnt.mean(x),\
                               "sterr":dumpcnt.sterr(x),\
                               "size":len(x)} for x in d[k]]))

pprint(gd)
gd = transpose(gd)
pprint(gd)
dumpcnt.vbarplot(gd)
