#!#!/usr/bin/env python3
"""Usage:
    synmap-to-key.py (--infile <arg>) (--outfile <arg>)

    -i INFILE --infile INFILE    input file from synmap
    -o OUTFILE --outfile OUTPUT  output file to write to for preliminary gene key
"""

from docopt import docopt
import re
REGEX = re.compile(r"(contig|ctg|scaffold|UNKNOWN)")

def process_lines(input):
    simple_lines = []
    with open(input) as i:
        for line in i:
            if line.startswith("#"):
                continue
            query_seqid,\
            query_meta,\
            query_start,\
            query_end,\
            subj_seqid,\
            subj_meta,\
            subj_start,\
            subj_end,\
            evalue,\
            logvalue = line.strip().split('\t')

            query,\
            query_chr,\
            query_st,\
            query_end,\
            per_id = query_meta.split('||')

            subj,\
            subj_chr,\
            subj_st,\
            subj_end,\
            per_id = subj_meta.split('||')

            parsed_query_chr = parseChrName(query_chr)
            parsed_subj_chr = parseChrName(subj_chr)
            if parsed_query_chr != parsed_subj_chr:
                m = re.match(REGEX, parsed_query_chr)
                p = re.match(REGEX, parsed_subj_chr)
                if m is None and p is None:
                    continue

            simple_lines.append((query_chr,
                                 query_st,
                                 query_end,
                                 query,
                                 subj_chr,
                                 subj_st,
                                 subj_end,
                                 subj,
                                 "synmap"))
    return simple_lines


def parseChrName(chromosome):
    chromosome = re.sub('chromosome', '', chromosome)
    chromosome = re.sub('chr', '', chromosome)
    m = re.match(r"0(\w+)", chromosome)
    if m:
        return m.group(1)
    else:
        return chromosome

def main():
    arguments = docopt(__doc__)
    formatted_lines = process_lines(arguments['--infile'])

    with open(arguments['--outfile'], "w") as o:
        for line in formatted_lines:
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                .format(line[0],
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
