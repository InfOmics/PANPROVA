import sys

linei = 0
for line in open(sys.argv[1]):
	if line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T':
		print(linei, line)
	linei += 1

