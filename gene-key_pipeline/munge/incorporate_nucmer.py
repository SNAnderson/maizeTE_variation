#!/usr/bin/env python3
"""Usage:
    incorporate-nucmer.py (--key <arg>) (--last <arg>) (--mummer_query <arg>) (--mummer_subj <arg>) (--output <arg>) (--query <arg>) (--subject <arg>) [--config <arg>]

    -k KEY --key KEY                   the preliminary gene key you're adding onto
    -l LAST --last LAST
    -q QUERY --query QUERY             the query genotype you are processing
    -s SUBJ --subject SUBJ             the subject genotype you are processing
    -m MUMMER --mummer_query MUMMER    parsed nucmer output to read from
    -p MUMMER --mummer_subj MUMMER     parsed nucmer
    -o OUTPUT --output OUTPUT          output file name
    -c CONFIG --config CONFIG          YAML config file [default: config/config.yaml]

    Used to incorporate genes paired through mummer alignment into gene key
"""
import sys,yaml,os,re
import locuspocus as lp
from collections import defaultdict
from docopt import docopt
# System call to bedtools intersect,
# get intersecting genes make sure coverage is at least 20
# also make sure order <4

def format_chrom(chromosome):
    return chromosome.strip('chr')

def isSyntenic(geneKey, queryOrder, subjOrder, queryGene, subjGene, queryGenome, subjGenome):
    if queryGene in geneKey:
        if subjGene in geneKey[queryGene]:
            return 1

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
            if abs(subjGenome[homolog_lookup].feature_order - subj_order) <= 3:
                downstream_pass = 1
                break
            else:
                continue

        if downstream_pass == 1:
            break

    while i < query_order + 10:
        i = i + 1
        lookup_gene = queryOrder[str(queryGenome[queryGene].chrom)][i]

        if not lookup_gene:
            continue
        if lookup_gene not in geneKey:
            continue

        for homolog_lookup in geneKey[lookup_gene]:
            if abs(subjGenome[homolog_lookup].feature_order - subj_order) <= 3:
                upstream_pass = 1
                break
            else:
                continue

        if upstream_pass == 1:
            break

    if upstream_pass == 1 or downstream_pass == 1:
        return 1
    else:
        return 0

def storeKey(geneKeyFile, query, subj, phb47key=None):
    '''Store the preliminary gene as a dict'''

    gene_key = defaultdict(lambda: defaultdict(list))
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
            gene_key[queryGene][subjGene] = line.strip().split('\t')

    return gene_key

def stripGeneName(gene, phb47_conversion=None):
    regex = re.compile(r"_(T|P)\d+$")
    gene=regex.sub("", gene).lstrip()
    if gene[0:5] == "ZePHB":
        gene = ".".join([gene, "g"])
        #gene = phb47_conversion[gene]

    return gene

def parseMummer(gene_key,last_hits,query_order,subj_order,query,subj,queryGenome,subjGenome, phb47key=None):
    a = defaultdict(list)
    b = defaultdict(list)

    with open('query.out') as f:
        for line in f:
            query_chr,\
            query_st,\
            query_sp,\
            query,\
            mmr_chr,\
            mmr_st,\
            mmr_sp,\
            coord = line.strip().split()
            query = stripGeneName(query,  phb47_conversion=phb47key)
            a[query].append(coord)

    with open('subj.out') as f:
        for line in f:
            subj_chr,\
            subj_st,\
            subj_sp,\
            subj,\
            mmr_chr,\
            mmr_st,\
            mmr_sp,\
            coord = line.strip().split()
            subj = stripGeneName(subj, phb47_conversion=phb47key)

            key = format_chrom(mmr_chr) + "-" + mmr_st + ";" + mmr_sp
            b[key].append(subj)

    for query_gene in a:
        for subj_coord in a[query_gene]:
            for subj_gene in b[subj_coord]:
                if format_chrom(queryGenome[query_gene].chrom) == format_chrom(subjGenome[subj_gene].chrom):
                    syntenyResult = isSyntenic(gene_key,query_order,subj_order,query_gene,subj_gene,queryGenome,subjGenome)
                    if syntenyResult == 1:
                        if subj_gene in gene_key[query_gene]:
                            gene_key[query_gene][subj_gene][8] = ",".join([gene_key[query_gene][subj_gene][8], "mummer"])
                        else:
                            if query_gene in last_hits:
                                if subj_gene in last_hits[query_gene]:
                                    gene_key[query_gene][subj_gene] = [queryGenome[query_gene].chrom, queryGenome[query_gene].start, queryGenome[query_gene].end, query_gene, subjGenome[subj_gene].chrom, subjGenome[subj_gene].start, subjGenome[subj_gene].end, subj_gene, 'mummer']

    return gene_key

def translatePHB47():
    import gffutils
    lookup = defaultdict(str)
    db = gffutils.FeatureDB('/home/maize/shared/databases/genomes/Zea_mays/PHB47/Zmaysvar.PHB47v1.1.gene_exons.gff3.db', keep_order=True)
    for gene in db.features_of_type('gene', order_by='start'):
        for i in db.children(gene, featuretype='mRNA', order_by='start'):
            lookup[i['Name'][0]] = gene['Name'][0]
    return lookup

def storeOrder(genome):
    '''Helper function to get gene order by chromosome'''
    feats = genome.get_feature_list(feature='gene')
    gene_order=defaultdict(lambda: defaultdict(dict))
    for feat in feats:
        gene_order[genome[feat].chrom][genome[feat].feature_order] = feat
    return gene_order

def storeLastHits(lastfile, phb47key=None):
    '''store hits from the last alignments'''
    h = defaultdict(list)
    #phb47_test = lp.RefLoci("PHB47", basedir="/home/maize/shared/databases/minus80")
    with open(lastfile) as l:
        for line in l:
            if line.startswith('#'):
                continue
            queryid,\
            subjectid,\
            peridentity,\
            alignmentlength,\
            mismatches,\
            gapopens,\
            qstart,\
            qend,\
            sstart,\
            send,\
            evalue,\
            bitscore = line.strip().split('\t')
            #translated_subject = phb47key[subjectid]
            queryid = stripGeneName(queryid)
            h[queryid].append(subjectid)

    return h

def main():
    args=docopt(__doc__)
    queryGenome = lp.RefLoci(args['--query'].upper(), basedir="/home/maize/shared/databases/minus80")
    subjGenome = lp.RefLoci(args['--subject'].upper(), basedir="/home/maize/shared/databases/minus80")
    with open (args['--config']) as c:
        config = yaml.safe_load(c)
        query_bed = config['samples'][args['--query'].lower()]['rep_bed']
        subj_bed = config['samples'][args['--subject'].lower()]['rep_bed']

    if args['--subject'].upper() == "PHB47":
        phb47_key = translatePHB47()

    cmd = " ".join(["bedtools intersect -wa -wb -f 0.80 -a", query_bed, "-b", args['--mummer_query'], "> query.out"])
    cmd2 = " ".join(["bedtools intersect -wa -wb -f 0.80 -a", subj_bed, "-b", args['--mummer_subj'], "> subj.out"])
    #os.system(cmd)
    #os.system(cmd2)
    #sys.exit()
    print("reading in last hits...")
    last_hits = storeLastHits(args['--last'], phb47key=None)
    print("storing query order...")
    query_order = storeOrder(queryGenome)
    print("storing subj order...")
    subj_order = storeOrder(subjGenome)
    print("storing gene key...")
    gene_key = storeKey(args['--key'], args['--query'], args['--subject'], phb47key=None)
    print("getting mummer results...")

    results = parseMummer(gene_key,last_hits,query_order,subj_order,args['--query'],args['--subject'],queryGenome,subjGenome, phb47key=None)
    with open(args['--output'], "w") as o:
        for query in results:
            for subj in results[query]:
                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                    .format(results[query][subj][0],
                            results[query][subj][1],
                            results[query][subj][2],
                            results[query][subj][3],
                            results[query][subj][4],
                            results[query][subj][5],
                            results[query][subj][6],
                            results[query][subj][7],
                            results[query][subj][8]))
    #os.system("rm query.out")
    #os.system("rm subj.out")

if __name__ == '__main__':
    main()
