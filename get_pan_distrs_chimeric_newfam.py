"""Chimera-aware pangenomic distributions, new-family variant.

Same as get_pan_distrs_chimeric.py but with the opposite reading of a chimera:
instead of multi-family membership, every chimeric gene founds a new family of
its own.
Reads .gene_parents (vertical lineage) written by evolve.cpp
and .chimeras.csv and generate single-membership family files:
a genuinely chimeric gene (an acceptor that received material from at least
one other gene) leaves its vertical family and starts a fresh family that its
vertical descendants inherit; every other gene resolves to its vertical root
as usual, so the output is a partition (every gene in exactly one family).

Output files (oprefix taken from argv[3]):
{oprefix}.chimeras_newfam.gene_families
{oprefix}.chimeras_newfam.family_presence
{oprefix}.chimeras_newfam.pan_distribution

These mirror the formats of get_pan_distrs.py outputs. Run alongside
get_pan_distrs_chimeric.py to compare the two interpretations on one run.

Usage:
python3 get_pan_distrs_chimeric_newfam.py <prefix>.gene_parents <prefix>.chimeras.csv <output_prefix>
"""

import sys
import csv


# args
if len(sys.argv) != 4:
    sys.stderr.write(
        "usage: get_pan_distrs_chimeric_newfam.py "
        "gene_parents_file chimeras_csv_file output_prefix\n"
    )
    sys.exit(1)

gene_parents_file = sys.argv[1]
chimeras_csv_file = sys.argv[2]
oprefix = sys.argv[3]


# ---------------------------------------------------------------------------
# parse .gene_parents
# format: space-separated 4 ints per line:
#   child_genome_id child_gene_id parent_genome_id parent_gene_id
# the (-1, -1) sentinel marks a family root (genome-0 ancestral genes
# and HGT-acquired genes both use it).
# ---------------------------------------------------------------------------

gene_parents = {}
gene_list = set()
genome_list = set()

with open(gene_parents_file, "r") as f:
    for line in f:
        cc = line.strip().split(" ")
        if len(cc) < 4:
            continue
        ig = (int(cc[0]), int(cc[1]))
        og = (int(cc[2]), int(cc[3]))
        gene_parents[ig] = og
        gene_list.add(ig)
        gene_list.add(og)
        genome_list.add(int(cc[0]))
        genome_list.add(int(cc[2]))

genome_list.discard(-1)
gene_list.discard((-1, -1))


# ---------------------------------------------------------------------------
# parse .chimeras.csv into the set of chimera founders.
# a founder is an acceptor gene that received material from at least one
# other gene (donor != acceptor). self-donor rows (donor == acceptor,
# e.g. INTRA sub-dup, INTRA inversion) are ignored: an intra-gene event
# merges no second lineage, so -- as in get_pan_distrs_chimeric.py, whose
# BFS skips self-donor edges -- it does not change the gene's family.
# ---------------------------------------------------------------------------

donors_of = {}

with open(chimeras_csv_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        acc = (int(row["acceptor_genome_id"]), int(row["acceptor_gene_id"]))
        don = (int(row["donor_genome_id"]), int(row["donor_gene_id"]))
        donors_of.setdefault(acc, []).append(don)

founders = set()
for acc, dons in donors_of.items():
    if any(d != acc for d in dons):
        founders.add(acc)


# ---------------------------------------------------------------------------
# single-family membership via vertical walk. returns the one family label
# of g: walk up the vertical-parent chain and stop at the first node that is
# a chimera founder (which cuts the chain and starts a new family) or a
# (-1, -1) root. founders are tested before the vertical parent is followed,
# so a chimeric gene founds its own family even though it still records a
# vertical parent. every node on the walked path shares that family, so the
# whole path is memoised in one shot. memo caches the result per gene.
# ---------------------------------------------------------------------------

def get_family(g, memo):
    path = []
    node = g
    seen = set()
    while True:
        if node in memo:
            fam = memo[node]
            break
        if node in seen:  # safety net against unexpected cycles
            fam = f"family_{node[0]}_{node[1]}"
            path.append(node)
            break
        seen.add(node)
        if node in founders:
            fam = f"family_chimera_{node[0]}_{node[1]}"
            path.append(node)
            break
        p = gene_parents.get(node)
        if p == (-1, -1) or p is None:
            fam = f"family_{node[0]}_{node[1]}"
            path.append(node)
            break
        path.append(node)
        node = p
    for n in path:
        memo[n] = fam
    return fam


family_to_genes = {}
family_to_genomes = {}
memo = {}

for g in gene_list:
    f = get_family(g, memo)
    family_to_genes.setdefault(f, set()).add(g)
    family_to_genomes.setdefault(f, set()).add(g[0])


# ---------------------------------------------------------------------------
# output: .chimeras_newfam.gene_families  (same format as .gene_families)
# ---------------------------------------------------------------------------

family_list = sorted(family_to_genes.keys())
genomes_sorted = sorted(genome_list)

with open(f"{oprefix}.chimeras_newfam.gene_families", "w") as off:
    for f in family_list:
        off.write(f + " ")
        for g in sorted(family_to_genes[f]):
            off.write(f"({str(g[0])},{str(g[1])}) ")
        off.write("\n")


# ---------------------------------------------------------------------------
# output: .chimeras_newfam.family_presence  (same format as .family_presence)
# ---------------------------------------------------------------------------

with open(f"{oprefix}.chimeras_newfam.family_presence", "w") as off:
    off.write("# " + " genome_".join(str(i) for i in genomes_sorted) + "\n")
    for f in family_list:
        off.write(f + " ")
        for g in genomes_sorted:
            off.write("x " if g in family_to_genomes[f] else "- ")
        off.write("\n")


# ---------------------------------------------------------------------------
# output: .chimeras_newfam.pan_distribution  (same format as .pan_distribution)
# ---------------------------------------------------------------------------

pand = {}
for f, gs in family_to_genomes.items():
    c = len(gs)
    pand[c] = pand.get(c, 0) + 1

with open(f"{oprefix}.chimeras_newfam.pan_distribution", "w") as off:
    off.write("#genomes #families\n")
    if pand:
        for c in range(1, max(pand.keys()) + 1):
            off.write(f"{str(c)} {str(pand.get(c, 0))}\n")
