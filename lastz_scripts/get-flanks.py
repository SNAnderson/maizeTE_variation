#!/usr/bin/env python3
import argparse
import sys
import locuspocus as lp

"""
ABOUT: This script is used to get the anchor genes for every TE throughout the genome.
RUNLIKE: % python munge/get-flanks.py --query-genotype b73 --subject-genotype w22 --outfile foo.out --dbpath /home/maize/shared/databases/minus80/
"""

def formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, rightHomolog, subject_loci):
    '''Format strings for the output format'''

    if "ctg" in query_locus.chrom or "scaffold" in query_locus.chrom or query_locus.chrom == "10000001" or query_locus.chrom == "UNKNOWN":
        if leftFlank is None:
            leftFlank = query_loci[rightFlank.chrom]
        if rightFlank is None:
            rightFlank = query_loci[leftFlank.chrom]
        homologPrint = "\t".join([str("NA") + "-" + str("NA"), str("NA")])
    else:
        if rightFlank.feature_type == "contig" or rightFlank.feature_type == "scaffold" or rightFlank.name == "UNKNOWN" or rightFlank.name =="10000001":
            homologPrint = "\t".join([str(subject_loci[leftHomolog].name) + "-" + str("NA"), str(subject_loci[leftHomolog].chrom) + ":" + str(subject_loci[leftHomolog].start) + "-" + str("NA")])
        elif leftFlank.feature_type == "contig" or leftFlank.feature_type == "scaffold" or leftFlank.name == "UNKNOWN" or leftFlank.name =="10000001":
            homologPrint = "\t".join([str("NA") + "-" + str(subject_loci[rightHomolog].name), str("NA") + ":" + str("NA") + "-" + str(subject_loci[rightHomolog].end)])
        else:
            homologPrint = "\t".join([str(subject_loci[leftHomolog.name].name) + "-" + str(subject_loci[rightHomolog.name].name), str(subject_loci[leftHomolog.name].chrom) + ":" + str(subject_loci[leftHomolog.name].start) + "-" + str(subject_loci[rightHomolog.name].end)])

    queryPrint = "\t".join([str(query_locus.name),
                            str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end)])
    flankPrint = "\t".join([str(leftFlank.name) + "-" + str(rightFlank.name),
                            str(leftFlank.chrom) + ":" + str(leftFlank.start) + "-" + str(rightFlank.end)])

    return "\t".join([queryPrint, flankPrint, homologPrint])

def getArgs():
    ''' Get arguments from command line '''
    parser = argparse.ArgumentParser(
    description='Script finds flanking loci for TE project')
    # Add arguments
    parser.add_argument('-q', '--query-genotype', required=True)
    parser.add_argument('-s', '--subject-genotype', required=True)
    parser.add_argument('-o', '--outfile', required=True)
    parser.add_argument('-p', '--dbpath', required=True)
    # Array for all arguments passed to script.
    args = parser.parse_args()
    # Assign args to variables
    query_genotype=args.query_genotype
    subject_genotype=args.subject_genotype
    outfile=args.outfile
    dbpath=args.dbpath
    # Return all variable values
    return query_genotype, subject_genotype, outfile, dbpath

def getFlanks(query, subject, dbpath):
    '''For each TE in the reference genome, get the closest flanking genes to serve as a window'''

    # This will load the database housing all of the relevant information
    query_loci = lp.RefLoci(query.upper(), basedir=dbpath)
    subject_loci = lp.RefLoci(subject.upper(), basedir=dbpath)
    printStatements = []

    # Now we loop over all of the loci
    for query_locus in query_loci:
        # We only care about TEs in this case
        if query_locus.feature_type != "transposable_element":
            continue
        # Get the closest flanking genes
        flanks = query_loci.flanking_loci(query_locus, window_size=100000000, feature="gene")
        ####################################################################################################
        ## If this gene is at the end or beginning of the chromosome we will have to treat it differently ##
        ####################################################################################################
        if len(flanks) != 2:
            # if there aren't any flanks then just print coordinates and NA -- should only happen for very small contigs
            if not flanks:
                printStatements.append("\t".join([str(query_locus.name), str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end), str("NA"), str("NA"), str("NA"), str("NA")]))
                continue

            elif flanks[0].feature_order > query_locus.feature_order:
                # We have at least one flank, this is if it's the right...
                rightFlank = flanks[0]
                if "ctg" in query_locus.chrom or "scaffold" in query_locus.chrom or query_locus.chrom == "    10000001" or query_locus.chrom == "UNKNOWN":
                    continue
                # If left flank is none then we make it equal to the chromosome name
                leftFlank = query_loci[rightFlank.chrom]
                leftHomolog = query_loci[rightFlank.chrom]

                # Try to get the homologs for this flanking gene, if there aren't any- move to the next outer flank
                rightFlank, rightHomolog = getHomologs(rightFlank, query_locus, query_loci, subject_loci, subject.upper(), "right")

                # No homologs found, probably at end of chrom. Add to results and continue
                if rightFlank is None:
                    printStatements.append("\t".join([str(query_locus.name), str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end), str("NA"), str("NA"), str("NA"), str("NA")]))
                    continue

                printStatements.append(formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, subject_loci[rightHomolog], subject_loci))
                continue

            else:
                # We have at least one flank, this if it's the left...
                leftFlank = flanks[0]

                if "ctg" in query_locus.chrom or "scaffold" in query_locus.chrom or query_locus.chrom == "10000001" or query_locus.chrom == "UNKNOWN":
                    continue
                # If right flank is none then we make it equal to the chromosome name
                rightFlank = query_loci[leftFlank.chrom]
                rightHomolog = query_loci[leftFlank.chrom]

                leftFlank, leftHomolog = getHomologs(leftFlank, query_locus, query_loci, subject_loci, subject.upper(), "left")

                # If you've moved off the edge of the chromsome ( should only happen for contigs )
                if leftFlank is None:
                    printStatements.append("\t".join([str(query_locus.name), str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end), str("NA"), str("NA"), str("NA"), str("NA")]))
                    continue
                printStatements.append(formatPrint(query_locus, leftFlank, subject_loci[leftHomolog], rightFlank, rightHomolog, subject_loci))
                continue

        #############################
        ## If we have both flanks ##
        #############################
        # Separate the flanking genes
        leftFlank, rightFlank = flanks

        ## If this a gene on a conting or unmapped chromosome, it won't be a good anchor
        if "ctg" in query_locus.chrom or "scaffold" in query_locus.chrom or query_locus.chrom == "10000001" or query_locus.chrom == "UNKNOWN":
            #printStatements.append(formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, rightHomolog, subject_loci))
            continue

        ###############################################
        ## Otherwise we should be on a 'normal' gene ##
        ###############################################
        # Flank may be None if we reached chromsome end when we searching for outer flank in previous block
        if leftFlank is None:
            leftHomolog = None
        else:
            # Otherwise, we get the homolog for this flank.
            leftFlank, leftHomolog = getHomologs(leftFlank, query_locus, query_loci, subject_loci, subject.upper(), "left")
        if rightFlank is None:
            rightHomolog = None
        else:
            rightFlank, rightHomolog = getHomologs(rightFlank, query_locus, query_loci, subject_loci, subject.upper(), "right")

        if leftHomolog is not None and rightHomolog is not None:
            # Similar to the flanks, we have to make sure the homologs are not the same gene
            if leftHomolog == rightHomolog:
                # Loop until a valid match is found
                while leftHomolog == rightHomolog:
                    # Get the next flanks and their homologs
                    leftFlank = getNextFlank(leftFlank, query_loci, query_locus, "left")
                    rightFlank = getNextFlank(rightFlank, query_loci, query_locus, "right")
                    if leftFlank is None or rightFlank is None:
                        break
                    leftFlank, leftHomolog = getHomologs(leftFlank, query_locus, query_loci, subject_loci, subject.upper(), "left")
                    rightFlank, rightHomolog = getHomologs(rightFlank, query_locus, query_loci, subject_loci, subject.upper(), "right")
                    if leftFlank is None or rightFlank is None or leftHomolog is None or rightHomolog is None:
                        break

                if leftHomolog is not None:
                    leftHomolog = subject_loci[leftHomolog]
                else:
                    leftFlank = query_loci[rightFlank.chrom]
                    leftHomolog = query_loci[rightFlank.chrom]
                if rightHomolog is not None:
                    rightHomolog = subject_loci[rightHomolog]
                else:
                    rightFlank = query_loci[leftFlank.chrom]
                    rightHomolog = query_loci[leftFlank.chrom]
                if leftFlank is None:
                    leftFlank = query_loci[rightFlank.chrom]
                if rightFlank is None:
                    rightFlank = query_loci[leftFlank.chrom]
                printStatements.append(formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, rightHomolog, subject_loci))
                continue

            # If the left homolog is downstream of the right homolog for these flanks, move to outer
            if subject_loci[leftHomolog].feature_order > subject_loci[rightHomolog].feature_order:
                # Loop until a valid match is found
                while subject_loci[leftHomolog].feature_order > subject_loci[rightHomolog].feature_order:
                    # Get the next flanks and their homologs
                    leftFlank = getNextFlank(leftFlank, query_loci, query_locus, "left")
                    rightFlank = getNextFlank(rightFlank, query_loci, query_locus, "right")
                    if leftFlank is None:
                        leftHomolog = None
                        break
                    if rightFlank is None:
                        rightHomolog = None
                        break
                    leftFlank, leftHomolog = getHomologs(leftFlank, query_locus, query_loci, subject_loci, subject.upper(), "left")
                    rightFlank, rightHomolog = getHomologs(rightFlank, query_locus, query_loci, subject_loci, subject.upper(), "right")
                    if leftFlank is None or rightFlank is None or leftHomolog is None or rightHomolog is None:
                        break
                if leftHomolog is not None:
                    leftHomolog = subject_loci[leftHomolog]
                else:
                    leftFlank = query_loci[rightFlank.chrom]
                    leftHomolog = query_loci[rightFlank.chrom]
                if rightHomolog is not None:
                    rightHomolog = subject_loci[rightHomolog]
                else:
                    rightFlank = query_loci[leftFlank.chrom]
                    rightHomolog = query_loci[leftFlank.chrom]
                printStatements.append(formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, rightHomolog, subject_loci))
                continue
        elif leftFlank is None:
            if rightFlank is not None:
                print(rightFlank)
                leftFlank = query_loci[rightFlank.chrom]
                leftHomolog = query_loci[rightFlank.chrom]
                rightFlank, rightHomolog = getHomologs(rightFlank, query_locus, query_loci, subject_loci, subject.upper(), "right")
                printStatements.append(formatPrint(query_locus, leftFlank, leftHomolog, rightFlank, subject_loci[rightHomolog], subject_loci))
                continue
            else:
                printStatements.append("\t".join([str(query_locus.name), str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end), str("NA"), str("NA"), str("NA"), str("NA")]))
                continue

        elif rightFlank is None:
            if leftFlank is not None:
                rightFlank = query_loci[leftFlank.chrom]
                rightHomolog = query_loci[leftFlank.chrom]
                leftFlank, leftHomolog = getHomologs(leftFlank, query_locus, query_loci, subject_loci, subject.upper(), "left")
                printStatements.append(formatPrint(query_locus, leftFlank, subject_loci[leftHomolog], rightFlank, rightHomolog, subject_loci))
                continue
            else:
                printStatements.append("\t".join([str(query_locus.name), str(query_locus.chrom) + ":" + str(query_locus.start) + "-" + str(query_locus.end), str("NA"), str("NA"), str("NA"), str("NA")]))
                continue

        printStatements.append(formatPrint(query_locus, leftFlank, subject_loci[leftHomolog], rightFlank, subject_loci[rightHomolog], subject_loci))
        print(leftFlank, rightFlank)
    return printStatements

def getHomologs(flank, queryLocus, queryGenome, subjGenome, subject, direction):
    '''Checks that both flanking genes in the query genotype have single homolog and that their homologs'''
    #if subject == "MO17":
    #    subject = "Mo17"
    # We require that a flank have a matched homolog in the opposite genotype
    rep_homolog = "NA"
    while rep_homolog == "NA":
        try:
            homologs = queryGenome.homologs(flank.name)
        except:
            homologs = ['.']
        if subject not in homologs:
            # If this gene didn't have any homologs, then we need to move to the next outer gene
            while subject not in homologs:
                # Get the next flank
                nextFlank = getNextFlank(flank, queryGenome, queryLocus, direction)
                if nextFlank is None:
                    return None, None
                # update homologs and flank
                homologs = queryGenome.homologs(nextFlank.name)
                flank = nextFlank

        rep_homolog = homologs[subject][0]
        # If there is more than one homolog, make sure we chose the "outer" most
        if len(homologs[subject]) > 1:
            for homolog in homologs[subject]:
                if direction == "left":
                    if int(subjGenome[homolog].feature_order) < int(subjGenome[rep_homolog].feature_order):
                        rep_homolog = homolog
                        order = int(subjGenome[homolog].feature_order)
                else:
                    if int(subjGenome[homolog].feature_order) > int(subjGenome[rep_homolog].feature_order):
                        rep_homolog = homolog
                        order = int(subjGenome[homolog].feature_order)

    return flank, rep_homolog

def getNextFlank(flank, queryGenome, queryLocus, direction):
    '''Get the next outer flank in the query genome'''
    passing = 0
    while passing != 1:
        if direction == "left":
            newFlank = queryGenome.upstream_loci(flank, window_size=100000000, feature="gene", locus_limit=1, within=1)
        else:
            newFlank = queryGenome.downstream_loci(flank, window_size=100000000, feature="gene", locus_limit=1, within=1)
        if not newFlank:
            return None

        # Avoid getting a new locus that's within then previous flank
        if queryLocus in queryGenome.loci_within(newFlank[0]):
            flank = newFlank[0]
            continue
        else:
            return newFlank[0]

def main():
    '''Main subroutine that makes calls to all relevant functions'''
    query, subject, outfile, dbpath = getArgs()
    printStatements = getFlanks(query, subject, dbpath)
    with open(outfile, 'w') as o:
        for printStatement in printStatements:
            o.write("%s\n" % str(printStatement))

if __name__ == '__main__':
    main()
