import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import statistics


i_evolve_log = "evolve.log"

df = pd.DataFrame()

print("gene counts")
gene_counts = list()
for line in open(i_evolve_log,'r'):
    if line.startswith("the new genome has a total of"):
        c = line.split(' ')[7]
        gene_counts.append(int(c))
print(min(gene_counts), max(gene_counts), statistics.mean(gene_counts), statistics.stdev(gene_counts))

df['gene_counts'] = gene_counts

print("genome lengths")
genome_lengths = list()
for line in open(i_evolve_log,'r'):
    if line.startswith("the length of the new genome is"):
        c = line.split(' ')[7]
        genome_lengths.append(int(c))
print(min(genome_lengths), max(genome_lengths), statistics.mean(genome_lengths), statistics.stdev(genome_lengths))

df['genome_lengths'] = genome_lengths

#plt.figure(figsize=(2, 4), dpi=80)
plt.figure(figsize=(6,6), dpi=80)
sns.boxplot(y='gene_counts', data=df)
sns.stripplot(y='gene_counts', data=df, color="orange", jitter=0.2, size=2.5)
plt.show()


#plt.figure(figsize=(2, 4), dpi=80)
plt.figure(figsize=(6,6), dpi=80)
sns.boxplot(y='genome_lengths', data=df)
sns.stripplot(y='genome_lengths', data=df, color="orange", jitter=0.2, size=2.5)
plt.show()



pand = dict()

i_pan_log = "pandistr.log"
dis = True
for line in open(i_pan_log,'r'):
    if dis :
        if line[0] == '=':
            dis = False
    else:
        cc = line.strip().split(' ')
        pand[ int(cc[0])] = int(cc[1])

print(pand)

pand_x = sorted(pand.keys())
pand_y = [ pand[x] for x in sorted(pand.keys()) ]


df = pd.DataFrame()
df['genomes'] = pand_x
df['genes'] = pand_y

plt.figure(figsize=(6,6), dpi=80)
#df.plot()
plt.plot(pand_x, pand_y)
plt.yscale('log')
plt.xlabel('N. of genomes', fontsize=10)
plt.ylabel('N. of genes', fontsize=10)
plt.show()
