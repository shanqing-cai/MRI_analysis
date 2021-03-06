#!/usr/bin/python

# Author: Shanqing Cai (shanqing.cai@gmail.com)
# Date: 12/20/2012


import os
import sys
import argparse
import tempfile
from subprocess import Popen, PIPE

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Rotate DWI gradient vectors (bvecs) to accomodate the eddy-current/motion correction parameters calculated by FSL eddy_corre.\nAssume:\n\t1.The input bvecs file is organized such that each line contains the three space-separated numbers, corresponding to the x, y and z components of a gradient.\n\t2. The ecclog file is in the format such as that generated by eddy_correct of FSL v4.1 or v5.0.\n\t3. FSL command avscale is available.")
    ap.add_argument("inBvecs", help="Input bvecs file")
    ap.add_argument("outBvecs", help="Output (rotated) bvecs file")
    ap.add_argument("ecclog", help="ecclog file generated by eddy_correct")
    ap.add_argument("-v", dest="bv", action="store_true", \
                    help="Verbose")

    if len(sys.argv) <= 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    inBvecs = args.inBvecs
    outBvecs = args.outBvecs
    ecclog = args.ecclog
    bv = args.bv

    # Check files
    if not os.path.isfile(inBvecs):
        raise Exception, "Cannot find input bvecs file: %s" % inBvecs
    if not os.path.isfile(ecclog):
        raise Exceptoin, "Cannot find ecclog file: %s" % ecclog

    # Check FSL commands
    (sout, serr) = Popen(["which", "avscale"], stdout=PIPE, stderr=PIPE)\
                   .communicate()
    if len(sout) == 0:
        raise Exception, "It appears that the required command avscale is unavailable."
    
    # Read input bvecs
    ibv = []

    bvf = open(inBvecs, "r")
    bvtxt = bvf.read().replace('\r', '\n').split('\n')
    bvf.close()

    for (i0, t_line) in enumerate(bvtxt):
        if len(t_line) == 0:
            continue

        if t_line.endswith(" "):
            t_line = t_line[:-1]

        if t_line.count(' ') != 2:
            raise Exception, "Unrecognized format in line %d of input bvecs file %s" (i0 + 1, inBvecs)
        tbv = []
        for word in t_line.split(' '):
            tbv.append(float(word))
        ibv.append(tbv)

    ngd = len(ibv)

    if bv:
        print("Found %d gradient directions." % ngd)

    # Read ecclog
    elf = open(ecclog, "r")
    eltxt = elf.read().replace('\r', '\n').split('\n')
    elf.close()
    
    # Rotation calculations
    obv = []
    tmpname = tempfile.mktemp()
    [foopath, basename] = os.path.split(tmpname)
    tmpfn = outBvecs + "_" + basename

    mats = []
    for i0 in range(ngd):
        t_obv = [0] * 3

        ln0 = 3 + i0 * 8
        ln1 = 6 + i0 * 8

        tmpf = open(tmpfn, "w")
        for i1 in range(ln0, ln1 + 1):
            t_line = eltxt[i1]
            if t_line.count(' ') != 4:
                raise Exception, "Unrecognized format in line #%d of input ecclog file: %s" % (i1 + 1, ecclog)
            tmpf.write(t_line + "\n")
        tmpf.close()
        
        (sout, serr) = Popen(["avscale", tmpfn], stdout=PIPE, stderr=PIPE)\
                       .communicate()
        avtxt = sout.replace('\r', '\n').split('\n')
        if avtxt[0] != "Rotation & Translation Matrix:" \
           or len(avtxt) < 20:
            raise Exception, "Unrecognized output format of avscale"
        m1 = avtxt[1].split(' ')
        m2 = avtxt[2].split(' ')
        m3 = avtxt[3].split(' ')
        m11 = float(m1[0])
        m12 = float(m1[1])
        m13 = float(m1[2])
        m21 = float(m2[0])
        m22 = float(m2[1])
        m23 = float(m2[2])
        m31 = float(m3[0])
        m32 = float(m3[1])
        m33 = float(m3[2])
        
        if bv:
            print("Frame #%d: rotation matrix = " % (i0 + 1))
            print("\t[%f, %f, %f;" % (m11, m12, m13))
            print("\t %f, %f, %f;" % (m21, m22, m23))
            print("\t %f, %f, %f]" % (m31, m32, m33))

        t_obv[0] = m11 * ibv[i0][0] + m12 * ibv[i0][1] + m13 * ibv[i0][2]
        t_obv[1] = m21 * ibv[i0][0] + m22 * ibv[i0][1] + m23 * ibv[i0][2]
        t_obv[2] = m31 * ibv[i0][0] + m32 * ibv[i0][1] + m33 * ibv[i0][2]
        obv.append(t_obv)

        if bv:
            print("Rotation result:")
            print("[%f, %f, %f] --> [%f, %f, %f]\n" % \
                  (ibv[i0][0], ibv[i0][1], ibv[i0][2], 
                   t_obv[0], t_obv[1], t_obv[2]))   
        
    # Write to output file
    outf = open(outBvecs, "w")
    for t_vec in obv:
        outf.write("%f %f %f\n" % (t_vec[0], t_vec[1], t_vec[2]))
    outf.close()
    if not os.path.isfile(outBvecs):
        raise Exception, "It appears that writing to output bvecs file %s has failed" % outBvecs

    if bv:
        print("Output written to file %s" % outBvecs)
        
    # Clean up
    os.system("rm -f %s" % tmpfn)
        
            
    
