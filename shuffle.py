#!/usr/bin/python

import sys
import random
import os
from math import log

if __name__ == "__main__":
    fnames = sys.argv[1:]
    d = os.path.commonprefix(fnames)
    random.shuffle(fnames)

    zpad = int(log(len(fnames), 10))

    dc = os.path.join(d, "count")
    os.mkdir(dc)
    for i in range(len(fnames)):
        if (i>0):
            fn = (zpad - int(log(i, 10))) * "0" + str(i)
        else:
            fn = "0" * (zpad + 1)
        lname = os.path.join(dc, fn + "." + fnames[i].split(".")[-1])
        os.symlink(fnames[i], lname)
        print("ln -s \"%s\" \"%s\"" %("../" + os.path.basename(fnames[i]),\
                os.path.basename(lname)))

