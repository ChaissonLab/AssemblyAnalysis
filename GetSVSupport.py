#!/usr/bin/env python


import sys
for line in sys.stdin:
    vals=line.split()
    svLen = int(vals[2]) - int(vals[1])
    
    support= vals[-1]
    nSup = 0
    if (support == "."):
        nSup = 0
    else:
        supVals = [int(i) for i in support.split(",") ]
        for supVal in supVals:
            if 0.5 < supVal / svLen < 2:
                nSup+=1
    vals = vals[0:5] + [str(nSup)]
    sys.stdout.write("\t".join(vals) + "\n")

        
    
