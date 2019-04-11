f = open("dump","r")
w = open("pair", "w")
for line in f:
	s=line.split()
	w.write(s[0]+"\n")
