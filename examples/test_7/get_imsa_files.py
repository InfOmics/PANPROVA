import sys

ifamilies = sys.argv[1]
igenes = sys.argv[2]
oprefix = sys.argv[3]


for line in open(ifamilies,'r'):
    cc = line.strip().split(' ')
    if len(cc) > 2:
        print(cc[0])
        ofile = open(oprefix+cc[0]+".fna",'w')

        for c in cc[1:]:
            ci = c.replace('(','').replace(')','').split(',')
            #print(ci)

            for gline in open(igenes, 'r'):
                gc = gline.strip().split(' ')
                gcc = gc[0].strip().split(':')
                if (gcc[0] == ci[0]) and (gcc[1] == ci[1]):
                    #print(ci,gline)
                    ofile.write('>' + gcc[0]+':'+gcc[1]+':'+gcc[2]+'\n')
                    ofile.write(gc[1]+'\n')

        ofile.flush()
        ofile.close()

        print('-'*40)