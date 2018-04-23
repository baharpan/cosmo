#!/usr/bin/env python

from collections import Counter
import sys

with open(sys.argv[1], 'r') as bloomfile:
    outfile = open(sys.argv[2], 'w')
    repeatfile = open("repeats.txt", 'w')
    kmer = ""
    reads = ""
    at_start = True
    
    for line in bloomfile:
	if "@" in line:
		continue
	
        x = line.split()
        if x[2] == "*" :
		continue
        if at_start:
            kmer = x[0]
	
            reads = [x[2]]
            at_start = False
        
        elif x[0] == kmer:
            reads.append(x[2])
        else:
            outfile.write(kmer + "\t")

            uniq = Counter(reads)
            read_nums = [i[0] for i in uniq.items()]
            read_counts = uniq.items()
            
            read_nums.sort(key=int)
            
            first = True
            for i in read_nums:
                if first:
                    outfile.write(i)
                    first = False
                else: 
                    outfile.write("," + i)
            outfile.write("\n")
            
            for i in read_counts:
                if i[1] > 1:
                    repeatfile.write(kmer + " ")
                    repeatfile.write(i[0] + "," + str(i[1]) + "\n")
            
            uniq.clear()
            kmer = x[0]
            del reads[:]
            del read_nums[:]
            reads = [x[2]]

    
    outfile.write(kmer + "\t")

    uniq = Counter(reads)
    read_nums = [i[0] for i in uniq.items()]
    read_counts = uniq.items()

    read_nums.sort(key=int)

    first = True
    for i in read_nums:
        if first:
            outfile.write(i)
            first = False
        else: 
            outfile.write("," + i)
    outfile.write("\n")
            
    
    for i in read_counts:
        if i[1] > 1:
            repeatfile.write(kmer + "\t")
            repeatfile.write(i[0] + "," + str(i[1]) + "\n")
    
    
    print("Extraction complete.\n")

repeatfile.close()
bloomfile.close()
outfile.close()
