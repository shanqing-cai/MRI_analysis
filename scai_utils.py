import os
import sys

def saydo(cmd, echo=True):
    if echo:
        print(cmd + '\n')

    os.system(cmd)

def qsubmit(cmd, queue, jobname):
    os.system('ezsub -c "%s" -q %s -n %s'%(cmd, queue, jobname))

def remove_empty_strings(lst):
    while lst.count('') > 0:
        lst.remove('')
    return lst

def read_ctab(ctabfn):
    ctabf = open(ctabfn, 'r')
    ctxt = ctabf.read().split('\n')
    ctabf.close()
    
    roi_nums = []
    roi_names = []
    for clin in ctxt:
        if len(clin) == 0:
            continue
        clin = clin.split(' ')
        while clin.count('') > 0:
            clin.remove('')
            
        if clin[1] == 'unknown' or clin[1] == 'bankssts' \
           or clin[1] == "corpuscallosum" or clin[1] == "Unknown" \
           or clin[1].startswith("None"):
            continue
        
        roi_nums.append(int(clin[0]))
        roi_names.append(clin[1])

    return (roi_nums, roi_names)

def check_file(fn):
    if not os.path.isfile(fn):
        raise Exception, "Missing file: %s"%fn

def check_dir(dn, bCreate=False):
    if not os.path.isdir(dn):
        if not bCreate:
            raise Exception, "Missing directory: %s"%dn
        else:
            os.system("mkdir %s"%dn)
            if not os.path.isdir(dn):
                raise Exception, "Failed to create directory: %s"%dn
            else:
                print("INFO: Created directory: %s"%dn)
