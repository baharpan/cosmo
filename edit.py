import sys
f = open(sys.argv[1],"r")
w = open(sys.argv[2], "w")
for line in f:
	s=line.split()
	w.write(s[0]+"\n")
