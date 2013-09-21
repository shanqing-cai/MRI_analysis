#!/usr/bin/python

import os
import sys
import numpy as np
import argparse

from scai_utils import *

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Swap rows and columns in an ASCII data file of numbers (e.g., diffusion-weighted MRI gradient orientation vectors")
    ap.add_argument("inFN", help="Input file name")
    ap.add_argument("outFN", help="Output file name")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    check_file(args.inFN)

    f = open(args.inFN, "rt")
    txt = remove_empty_strings(f.read().split("\n"))
    f.close()

    nrows = len(txt)
    
    a_ncols = np.zeros(nrows)
    for (i0, tline) in enumerate(txt):
        items = remove_empty_strings(tline.split(" "))
        a_ncols[i0] = len(items)

    assert(len(np.unique(a_ncols)) == 1)

    ncols = int(a_ncols[0])
    dat = np.zeros([nrows, ncols])

    for (i0, tline) in enumerate(txt):
        items = remove_empty_strings(tline.split(" "))
        
        for (i1, t_item) in enumerate(items):
            dat[i0, i1] = float(t_item)
    

    dat = np.transpose(dat)

    
    f = open(args.outFN, "wt")
    for i0 in range(ncols):
        for i1 in range(nrows):
            f.write("%.14f " % dat[i0, i1])
        f.write("\n")

    f.close()
    check_file(args.outFN)
