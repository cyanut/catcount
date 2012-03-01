#!/usr/bin/python

import xml.etree.ElementTree as etree
import sys
import os

import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import matplotlib.ticker as ticker
import numpy as np
from collections import OrderedDict
from math import sqrt
import copy

def getcounts(flist):
    countres={}
    for fname in flist:
        imgname = ""
        res = []
        tree = etree.parse(fname)
        root = tree.getroot()
        imgname = root.getchildren()[0].getchildren()[0].text
        imgname = os.path.realpath(os.path.join(os.path.dirname(fname),imgname))
        mdata = root.getchildren()[1]
        for i in range(len(mdata.getchildren()) - 1):
            res.append(-1)
        for mt in mdata.getchildren():
            if mt.tag == "Marker_Type":
                cl = mt.getchildren()
                n = int(cl[0].text)
                res[n - 1] = len(cl) - 1
        countres[imgname] = res
    return countres

def mean(l):
    return sum(l)/1.0/len(l)

def sterr(l):
    if len(l) == 1:
        return 0
    m = mean(l)
    return sqrt(sum([(i-m) * (i-m) for i in l])/(len(l) - 1))/sqrt(len(l))  

def poolgroup(animalres, groupdic):
    groupdata = OrderedDict()
    for gn in groupdic:
        groupdata[gn] = OrderedDict()
        tgd = OrderedDict()
        for an in groupdic[gn]:
            if an in animalres:
                for mn in animalres[an]:
                    if not mn in tgd:
                        tgd[mn] = [animalres[an][mn]]
                    else:
                        tgd[mn].append(animalres[an][mn])
        for mn in tgd:
            groupdata[gn][mn] = OrderedDict([("mean", mean(tgd[mn])),\
                    ("sterr", sterr(tgd[mn])), ("size", len(tgd[mn]))])
    return groupdata

def getgroupdic(fname):
    groupname = {}
    groupdic = OrderedDict()
    f = open(fname)
    for line in f:
        lt =line.strip()
        if len(lt) > 0 and lt[0] == "#":
            continue
        l = line.split("|")
        if len(l) >= 2:
            if l[0][0] == "G":
                groupname[l[0].strip()] = l[1].strip()
                groupdic[l[1].strip()] = []
            else:
                groupdic[groupname[l[1].strip()]].append(l[0].strip())
    f.close()
    return groupdic

def getcomments(fname):
    commentdic = {}
    f = open(fname)
    for line in f:
        l = line.split("|")
        if len(l) > 2:
            commentdic[l[0].strip()] = l[2].strip()
    f.close()
    return commentdic

def getanimalcode(imgpath):
    fname = os.path.basename(imgpath)
    animcode = re.findall("^(.+?[0-9]+)", fname)
    if animcode:
        return ''.join(x for x in animcode[0] if x != '-')
    else:
        return None

def poolanimal(countres, postfunc=""):
    if postfunc == "":
        def postfunc(m):
            return m
    animalres = {}
    finalanimalres = {}
    for imgpath in countres:
        n = getanimalcode(imgpath)
        if not n:
            continue
        if (not n in animalres):
            animalres[n] = copy.copy(countres[imgpath])
        else:
            la = len(animalres[n])
            lc = len(countres[imgpath])
            for i in range(min(lc, la)):
                animalres[n][i] += countres[imgpath][i]
            for i in range(lc-la):
                animalres[n].append(countres[imgpath][i + la])
    for n in animalres:
        finalanimalres[n] = postfunc(animalres[n])
    return finalanimalres


def count_arc(cellcount):
    return OrderedDict([\
            ("DAPI", cellcount[0]), \
            ("Arc", cellcount[1]), \
            ("GFP", cellcount[2]), \
            ("GFP/Arc", cellcount[3])])

def prob_arc(cellcount):
    '''input: a list of counts [DAPI, Arc, GFP, Arc/GFP]
    output: a ordereddict of [pGFP, pArc, pArc/GFP-, pArc/GFP+]
    '''
    return OrderedDict([\
            #("pGFP", cellcount[2]/1.0/cellcount[0]), \
            ("Total", cellcount[1]/1.0/cellcount[0]), \
            ("GFP-", (cellcount[1]-cellcount[3])/1.0/\
            (cellcount[0]-cellcount[2])), \
            ("!GFP+", cellcount[3]/1.0/cellcount[2])])

def dirty_prob_arc(cellcount):    
    return OrderedDict([\
            #("pGFP", cellcount[2]/1.0/cellcount[0]), \
            ("Total", cellcount[1]/1.0/cellcount[0]), \
            ("GFP-", (cellcount[1]-cellcount[3])/1.0/\
            (cellcount[0]-cellcount[2])), \
            ("GFP+", cellcount[3]/1.0/cellcount[2]),\
            ("!pChance+", (cellcount[3] - 1.0*cellcount[1]/cellcount[0]*cellcount[2])/cellcount[0])])

def prob_arc_gfp(cellcount):
    return OrderedDict([\
            ("pGFP", cellcount[2]/1.0/cellcount[0]), \
            ("pArc", cellcount[1]/1.0/cellcount[0]), \
            ("pBoth", cellcount[3]/1.0/cellcount[0]), \
            ("pChance", cellcount[1]*cellcount[2]/1.0/cellcount[0]/cellcount[0])])

def getanimalnumber(imgpath, prefix="sl"):  
    fname = os.path.basename(imgpath)
    animnum = re.findall("^"+prefix+".*?([0-9]+)", fname)
    if animnum:
        return animnum[0]
    else:
        return None


def vbarplot(groupdata):
    xoffset = -0.2
    goffset = 0.4
    width = 0.95
    axiscolor = "black"
    barcolor = "#b9b9b9" 
    sigbarcolor = "#555555"
    errcolor = "black"
    ylabel = "p(Arc+)"
    ticks = []
    xlabel = []
    grpliney = -0.055
    grplinew = 1.5
    grplinespan = 0.3
    grptxty = -0.10
    barlw = 1.5
    errlw = barlw
    axislw = barlw

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for s in ax.spines.values():
        s.set_linewidth(axislw)
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    #ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.005))
    for gn,gd in groupdata.items():
        mname = gd.keys()
        mmean = [m["mean"] for m in gd.values()]
        merr = [m["sterr"] for m in gd.values()]
        ind = np.linspace(xoffset, xoffset + len(mname) - 1, len(mname))
        ticks.extend(list(ind))
        for i in range(len(ind)):
            bc = barcolor
            if mname[i][0] == '!':
                mname[i] = mname[i][1:]
                bc = sigbarcolor
            ax.bar(ind[i], mmean[i], width, align="center", color=bc, aa=True, lw=barlw)
            ax.errorbar(ind[i], mmean[i], yerr=merr[i], ecolor=errcolor, elinewidth=errlw)
        xlabel.extend(mname)
        xoffset += len(mname) + goffset 
        
    ax.autoscale(False, None, True)

    grply = ax.transData.inverted().transform(ax.transAxes.transform((0, grpliney)))[1]
    grptxty = fig.transFigure.inverted().transform(ax.transAxes.transform((0, grptxty)))[1]

    nitem = 0
    for gn,gd in groupdata.items():
        l = Line2D((ticks[nitem] - grplinespan, ticks[nitem + len(gd) -1] + grplinespan), (grply, grply), lw=grplinew, color="black", clip_on=False)
        ax.add_line(l)
        #igrp_text = Text((ticks[nitem] + ticks[nitem + len(gd) - 1])/2, grptxty, text=gn, ha="center", axes=ax)
        #ax.add_artist(grp_text)
        
        gtx = (ticks[nitem] + ticks[nitem + len(gd) -1])/2
        gtx = fig.transFigure.inverted().transform(ax.transData.transform((gtx, 0)))[0]
        gname ="%s (n=%d)" %(gn, groupdata[gn][groupdata[gn].keys()[0]]["size"])
        plt.figtext(gtx, grptxty, gname, ha="center", size="large", weight="bold")
        nitem += len(gd)
    ax.set_xticks(ticks)
    ax.set_xticklabels(xlabel, size="large")
    ax.set_ylabel(ylabel, weight="bold", size="large")
    plt.show()

def print_table(animalres, dic, comments={}):
    first = True
    gn = ""
    for key in animalres:
        if first:
            print("Animal #\tCondition\t" + "\t".join(animalres[key].keys()) +\
                    "\tComment")
            first = False
        for gname in dic:
            if key in dic[gname]:
                gn = gname.strip()
                print(" \t".join((key, gn, \
                "\t".join(["%.4f" %i for i in animalres[key].values()]),\
                comments.get(key, ""))))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: "+sys.argv[0]+" groupfile prefix imgfiles...")
        quit()
    groupf = sys.argv[1]
    imgcounts = getcounts(sys.argv[2:])
    for key in imgcounts:
        print(key + "|" + "|".join([str(i) for i in imgcounts[key]]))
    c = poolanimal(imgcounts, postfunc=prob_arc)
    counts = poolanimal(imgcounts, postfunc=count_arc)
    allc = poolanimal(imgcounts, postfunc=prob_arc_gfp)
    #dc = poolanimal(imgcounts, postfunc=dirty_prob_arc)
    d = getgroupdic(groupf)
    cm = getcomments(groupf)
    print_table(c, d, cm)
    print_table(counts, d, cm)
    print_table(allc, d, cm)
    #print_table(dc, d, cm)
    g = poolgroup(c, d)
    vbarplot(g)
