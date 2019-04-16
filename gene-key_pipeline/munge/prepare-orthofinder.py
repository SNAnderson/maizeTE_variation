# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""Usage:
    prepare-orthofinder.py (--query <arg>) (--subj <arg>) (--querybed <arg>) (--subjbed <arg>) (--dir <arg>) [--config <arg>]

    --query QUERY            run mode (fasta / bed)
    --subj SUBJ              genotype to analyze
    --querybed QUERYBED      querybed file with representative features
    --subjbed SUBJBED        subjbed file with representative features
    --config CONFIG          YAML config file [default: config/config.yaml]
    --dir CONFIG             YAML config file [default: config/config.yaml]
"""

import yaml,sys,os,re
from docopt import docopt
from pyfaidx import Fasta
from natsort import natsorted

def stripName(gene):
    """get representative cds id from bed file"""
    lookup = {}
    m = re.search('(\w+)_[T|P](\d+)', gene)
    if m is not None:
        base = m.group(1)
        child = m.group(2)
    else:
         base = child = gene
    return base, child

def remakeProt(fasta, outfile, idfile, id):
    fasta_index = Fasta(fasta)
    lookup = {}
    for protein in fasta_index.keys():
        size = len(fasta_index[protein])
        gene, isoform = stripName(protein)
        if gene in lookup:
            if size > lookup[gene][0]:
                lookup[gene] = (size, isoform)
        else:
            lookup[gene] = (size, isoform)
    
    with open(outfile,"w") as f, open(idfile, "a") as q:
        for i, gene in enumerate(lookup):
            isoform = lookup[gene][1]
            name = "".join([gene, "_P", isoform])
            if name not in fasta_index:
                if gene == isoform:
                    name = gene
                else:
                    name = "".join([gene, "_T", isoform])
            q.write("{}_{}: {}\n".format(id, i, gene))
            f.write(">{}_{}\n".format(id,i))
            for line in fasta_index[name]:
                f.write("{}\n".format(str(line)))

def main():
    arguments = docopt(__doc__)
    print(arguments)
    with open(arguments['--config']) as c:
       config = yaml.safe_load(c)
       query_fasta = config['samples'][arguments['--query']]['prot']
       subj_fasta = config['samples'][arguments['--subj']]['prot']
    

    cwd = os.getcwd()
    full_dir = os.path.join(cwd, arguments['--dir'])
    if not os.path.exists(full_dir):
        os.makedirs(full_dir)

    query_fasta_out = os.path.join(full_dir, "Species0.fa")
    subj_fasta_out = os.path.join(full_dir,  "Species1.fa")
    idfile = os.path.join(full_dir, 'SequenceIDs.txt')
    speciesfile = os.path.join(full_dir, 'SpeciesIDs.txt')
    remakeProt(query_fasta, query_fasta_out, idfile, '0')
    remakeProt(subj_fasta, subj_fasta_out, idfile, '1')
    with open(speciesfile, "w") as s:
        s.write("0: {}\n".format(query_fasta))
        s.write("1: {}".format(subj_fasta))

if __name__ == '__main__':
    main()
