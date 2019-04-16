#!/usr/bin/env python2
"""Usage:
    get-rep-cds.py (--genotype <arg>) (--mode fasta|bed) (--outfile <arg>) [--size <arg>] [--parent <arg>] [--child <arg>] [--extract_id <arg>] [--config <arg>]

    -q GENOTYPE --genotype GENOTYPE  query genotype to analyze
    -m MODE --mode mode              mode to run (fasta / bed)
    -o OUTFILE --outfile OUTFILE     output file name
    -p PARENT --parent PARENT        parent feature [default: gene]
    -d CHILD --child CHILD           child feature [default: mRNA]
    -f FEAT --extract_id FEAT        feature connecting CDS fasta to CDS id in GFF  [default: Parent]
    -c CONFIG --config CONFIG        YAML config file [default: config/config.yaml]
    -s SIZE --size SIZE              cutoff size [default: 10000]

    Used to classify TEs into various categories of being nearby, overlapping, or within genes.
"""

import yaml
import sys
import os
from docopt import docopt
from jcvi.formats.gff import *
from pyfaidx import Fasta
from natsort import natsorted
from collections import defaultdict

def stripName(name):
    m = re.search('(\w+)\.v1\.\d+\.CDS\.\d+', name)
    if m is not None:
        name = m.group(1)
    n = re.search('(\w+)\.v\d+\.\d+', name)
    if n is not None:
        name = n.group(1)
    o = re.search('(\w+)_CDS', name)
    if o is not None:
        name = o.group(1)
    parsed_name = name.split(':')
    if len(parsed_name) > 1:
        return parsed_name[1]
    else:
        return parsed_name[0]

def get_lengths(input, cutoff):
    lookup = {}
    contigs = Fasta(input)

    for contig in contigs.keys():
        if len(contigs[contig]) < cutoff:
            lookup[contig] = 1

    return lookup

def getReps(gff_index, fasta, excluded_contigs, outfile, extract_id, mode, parent, cftype):
    reps = defaultdict(lambda: defaultdict(dict))

    print("getting sizes...")
    for feat in gff_index.features_of_type(parent):
        children = [child for child in gff_index.children(feat, level = 1, featuretype=cftype)]
        for child in children:
            if feat.id not in reps:
                reps[feat.id] = (child.id, gff_index.children_bp(feature = child, child_featuretype = "CDS"))
            else:
                if gff_index.children_bp(feature = child, child_featuretype = "CDS") > reps[feat.id][1]:
                    reps[feat.id] = (child.id, gff_index.children_bp(feature = child, child_featuretype = "CDS"))

    if mode == "fasta":
        fasta_index = Fasta(fasta)

        with open(outfile,"w") as f:

            for rep_feature in reps:
                rep_id, size = reps[rep_feature]

                if gff_index[rep_id].chrom in excluded_contigs:
                    continue

                child_id = [c.id for c in gff_index.children(rep_feature, featuretype = "CDS")][0]
                try: 
                    fasta_id = gff_index[child_id][extract_id][0]
                except:
                    fasta_id = gff_index[gff_index[child_id]['Parent'][0]][extract_id][0]
                f.write(">{}\n".format(stripName(rep_feature)))
                for line in fasta_index[stripName(fasta_id)]:
                    f.write("{}\n".format(str(line)))

    elif mode == "bed":
        with open(outfile, "w") as o:
            for feat in gff_index.features_of_type("gene"):
                feat_rep = reps[feat.id][0]
                if gff_index[feat_rep].chrom in excluded_contigs:
                    continue

                o.write("{}\t{}\t{}\t{}\n".format(gff_index[feat_rep].chrom, gff_index[feat_rep].start-1, gff_index[feat_rep].end, stripName(feat.id)))

def main():
    arguments = docopt(__doc__)
    with open(arguments['--config']) as c:

       config = yaml.safe_load(c)
       gffile = config['samples'][arguments['--genotype']]['gff']
       fasta = config['samples'][arguments['--genotype']]['cds']
       genome = config['samples'][arguments['--genotype']]['fasta']

    gff = make_index(gffile)
    excluded_contigs = get_lengths(genome, int(arguments['--size']))
    getReps(gff, fasta, excluded_contigs, arguments['--outfile'], arguments['--extract_id'], arguments['--mode'], arguments['--parent'], arguments['--child'])
    pwd = os.getcwd()
    new_path = os.path.join(pwd, arguments['--outfile'])

    if arguments['--mode'] == "bed":
        config['samples'][arguments['--genotype']]['rep_bed'] = new_path
    elif arguments['--mode'] == "fasta":
        config['samples'][arguments['--genotype']]['rep_fasta'] = new_path
    with open(arguments['--config'], 'w') as y:
        yaml.dump(config,y)

if __name__ == '__main__':
    main()
