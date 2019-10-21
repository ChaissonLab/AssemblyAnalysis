#!/usr/bin/env python
import sys

title=["correct-het", "correct-hom", "missed-het", "missed-hom", "fp-het", "fp-hom", "het-to-hom-alt", "hom-alt-to-het"]
correct_het=0
correct_hom=1
missed_het=2
missed_hom=3
fp_het=4
fp_hom=5
het_to_hom_alt=6
hom_alt_to_het=7
n=8
tally = [[0]*n, [0]*n]



for line in sys.stdin:
    if line[0] == "#":
        continue
    vals=line.split()
    gts=[vals[-3].split(":")[0], vals[-2].split(":")[0]]
    sts=vals[-1].split(":")[0].replace("|","/")
#    print((gts,sts))
    if sts == "1/0":
        sts = "0/1"

    for i in [0,1]:
        if sts == "0/0":
            if gts[i] == "0/1":
                tally[i][fp_het]+=1
            elif gts[i] == "1/1":
                tally[i][fp_hom]+=1
        if sts == "0/1":
            if gts[i] == "0/1":
                tally[i][correct_het]+=1
            elif gts[i] == "1/1":
                tally[i][het_to_hom_alt]+=1
            elif gts[i] == "0/0":
                tally[i][missed_het]+=1
        if sts == "1/1":
            if gts[i] == "0/0":
                tally[i][missed_hom]+=1
            if gts[i] == "0/1":
                tally[i][hom_alt_to_het] +=1
            if gts[i] == "1/1":
                tally[i][correct_hom] +=1

print("\t".join(title))
print("\t".join([str(i) for  i in tally[0]]))
print("\t".join([str(i) for i in tally[1]]))
                
                
    
    
    
