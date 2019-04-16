#!/usr/bin/env python3
"""Usage:
    incorporate-orthofinder.py (--key <arg>) (--infile <arg>) [--outfile <arg>] (--query <arg>) (--subj <arg>)

    -k KEY --key KEY              gene-key file
    -i INFILE --infile INFILE     input file from synmap
    -q QUERY --query QUERY        query genotype
    -s SUBJ --subj SUBJ           subject genotype
    -o OUTFILE --outfile OUTFILE  output file to write to for preliminary gene key

RUNLIKE python3 munge/incorporate-orthofinder.py --infile b73-phb47/orthofinder/Orthogroups.csv --query B73 --subj PHB47 --key b73-phb47/b73-phb47_key.txt --outfile foo.txt

"""

from docopt import docopt
import locuspocus as lp
from collections import defaultdict
import sys,re

def format_chrom(chromosome):
    return chromosome.strip('chr')

def isSyntenic(geneKey, queryOrder, subjOrder, queryGene, subjGene, queryGenome, subjGenome):
    if queryGene in geneKey:
        print("query in key")
        if subjGene in geneKey[queryGene]:
            return 2

    query_order = queryGenome[queryGene].feature_order
    subj_order = subjGenome[subjGene].feature_order

    i = query_order
    upstream_pass = downstream_pass = 0

    while i > query_order - 10 and i > 0:
        i = i-1
        lookup_gene = queryOrder[queryGenome[queryGene].chrom][i]
        if not lookup_gene:
            continue

        if lookup_gene not in geneKey:
            continue
        for homolog_lookup in geneKey[lookup_gene]:
            if abs(subjGenome[homolog_lookup].feature_order - subj_order) <= 8:
                #print(lookup_gene,homolog_lookup)
                downstream_pass = 1
                break
            else:
                continue

        if downstream_pass == 1:
            break

    while i < query_order + 10:
        i = i + 1
        lookup_gene = queryOrder[queryGenome[queryGene].chrom][i]

        if not lookup_gene:
            continue
        if lookup_gene not in geneKey:
            continue

        for homolog_lookup in geneKey[lookup_gene]:
            if abs(subjGenome[homolog_lookup].feature_order - subj_order) <= 8:
                #print(lookup_gene,homolog_lookup)
                upstream_pass = 1
                break
            else:
                continue

        if upstream_pass == 1:
            break

    if upstream_pass == downstream_pass == 1:
        return 1
    else:
        return 0

def stripGeneName(hit, genotype, phb47_conversion=None):
    regex = re.compile(r"_(T|P)\d+$")
    hit=regex.sub("", hit).lstrip()
    if genotype.upper() == "PHB47":
        hit = phb47_conversion[hit]
    return hit

def parseOrthofinder(orthofile,gene_key,query_order,subj_order,query,subj,queryGenome,subjGenome,phb47key=None):
    final_result = defaultdict(list)
    seen = defaultdict()
    with open(orthofile) as o:
        next(o)
        regex = re.compile("^Zm")
        for line in o:
            fields = line.strip().split('\t')
            if len(fields) != 3:
                continue

            groups,\
            queryGenes,\
            subjGenes = fields

            queryHits = queryGenes.split(',')
            subjHits =  subjGenes.split(',')
            if re.search(regex, queryHits[0]) is None:
                continue
            if re.search(regex, queryHits[0]) is None:
                continue

            for queryHit in queryHits:
                for subjHit in subjHits:
                    queryHit = stripGeneName(queryHit, query, phb47_conversion=phb47key)
                    subjHit = stripGeneName(subjHit, subj, phb47_conversion=phb47key)

                    if "|".join([queryHit,subjHit]) in seen:
                        continue
                    else:
                        seen["|".join([queryHit,subjHit])] = 1
                    print("query:")
                    print(queryHit)
                    print(format_chrom(queryGenome[queryHit].chrom))
                    print("subj:")
                    print(subjHit)
                    print(format_chrom(subjGenome[subjHit].chrom))
                    if format_chrom(queryGenome[queryHit].chrom) == format_chrom(subjGenome[subjHit].chrom):
                        syntenyResult = isSyntenic(gene_key,query_order,subj_order,queryHit,subjHit,queryGenome,subjGenome)
                        print(syntenyResult)
                        if syntenyResult == 1:
                            final_result[queryHit].append((subjHit, 'orthofinder'))
                        elif syntenyResult == 2:
                            final_result[queryHit].append((subjHit, 'synmap,orthofinder'))
                        else:
                            continue
    return final_result

def storeKey(geneKeyFile, query, subj, phb47key=None):
    '''Store the preliminary gene as a dict'''

    gene_key = defaultdict(list)
    with open(geneKeyFile) as k:
        for line in k:
            queryChrom,\
            queryStart,\
            queryEnd,\
            queryGene,\
            subjChrom,\
            subjStart,\
            subjEnd,\
            subjGene,\
            evidence = line.strip().split('\t')

            queryGene = stripGeneName(queryGene, query, phb47_conversion=phb47key)
            subjGene = stripGeneName(subjGene, subj, phb47_conversion=phb47key)
            gene_key[queryGene].append(subjGene)
            gene_key[subjGene].append(queryGene)
        return gene_key

def storeOrder(genome):
    '''Helper function to get gene order by chromosome'''
    feats = genome.get_feature_list(feature='gene')
    gene_order=defaultdict(lambda: defaultdict(dict))
    for feat in feats:
        gene_order[genome[feat].chrom][genome[feat].feature_order] = feat
    return gene_order

def translatePHB47():
    import gffutils
    lookup = defaultdict(str)
    db = gffutils.FeatureDB('/home/maize/shared/databases/genomes/Zea_mays/PHB47/Zmaysvar.PHB47v1.1.gene_exons.gff3.db', keep_order=True)
    for gene in db.features_of_type('gene', order_by='start'):
        for i in db.children(gene, featuretype='mRNA', order_by='start'):
            lookup[i['Name'][0]] = gene['Name'][0]
    return lookup

def main():
    arguments = docopt(__doc__)

    queryGenome = lp.RefLoci(arguments['--query'].upper(), basedir="/home/maize/shared/databases/minus80")
    subjGenome = lp.RefLoci(arguments['--subj'].upper(), basedir="/home/maize/shared/databases/minus80")
    print("storing query order...")
    query_order = storeOrder(queryGenome)
    print("storing subj order...")
    subj_order = storeOrder(subjGenome)
    print("getting ortho results...")

    print("storing gene key...")
    if arguments['--subj'].upper() == "PHB47":
        phb47_key = translatePHB47()
        gene_key = storeKey(arguments['--key'], arguments['--query'], arguments['--subj'],phb47key=phb47_key)
        ortho_results = parseOrthofinder(arguments['--infile'],gene_key,query_order,subj_order,arguments['--query'],arguments['--subj'],queryGenome,subjGenome,phb47key=phb47_key)
    else:
        gene_key = storeKey(arguments['--key'], arguments['--query'], arguments['--subj'])
        ortho_results = parseOrthofinder(arguments['--infile'],gene_key,query_order,subj_order,arguments['--query'],arguments['--subj'],queryGenome,subjGenome)

    with open(arguments['--outfile'], "w") as o:
        feats = queryGenome.get_feature_list(feature='gene')
        for feat in feats:
            lines=[]
            if feat in gene_key and feat in ortho_results:
                for result in ortho_results[feat]:
                    hit, evidence = result
                    lines.append((queryGenome[feat].chrom,
                                  queryGenome[feat].start,
                                  queryGenome[feat].end,
                                  feat,
                                  subjGenome[hit].chrom,
                                  subjGenome[hit].start,
                                  subjGenome[hit].end,
                                  hit,
                                  evidence))
            if feat in gene_key and feat not in ortho_results:
                for hit in gene_key[feat]:
                    lines.append((queryGenome[feat].chrom,
                                  queryGenome[feat].start,
                                  queryGenome[feat].end,
                                  feat,
                                  subjGenome[hit].chrom,
                                  subjGenome[hit].start,
                                  subjGenome[hit].end,
                                  hit,
                                  str("synmap")))
            if feat not in gene_key and feat in ortho_results:
                for result in ortho_results[feat]:
                    hit, evidence = result
                    lines.append((queryGenome[feat].chrom,
                                  queryGenome[feat].start,
                                  queryGenome[feat].end,
                                  feat,
                                  subjGenome[hit].chrom,
                                  subjGenome[hit].start,
                                  subjGenome[hit].end,
                                  hit,
                                  evidence))
            for line in lines:
                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    line[0],
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                    line[5],
                    line[6],
                    line[7],
                    line[8]))

if __name__ == '__main__':
    main()
