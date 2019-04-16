#!/usr/bin/env python3
""" This script parses lav outputs and determines the coverage of all the TEs in the window and determine what TEs in the opposite
    genotype hit to.
    See https://www.bx.psu.edu/miller_lab/dist/lav_format.html for info on the LAV format.
ABOUT: This script will parse the lav formatted lastz outputs to determine the coverage for each element in the window.
It will also report whether it aligned to any other elements in the opposite genotype and the bestHit based solely on coverage.
"""

import sys,re,glob,os
from collections import defaultdict

def parse_d_stanza(l):
    'This basically just skips the d stanza, which is the scoring matrix used'

    line = l.readline()
    while line:
        line = l.readline()
        if line.startswith("}"):
            break
        else:
            continue

def parse_s_stanza(l):
    'The S-stanza contains info on the filenames, which are the flanking genes'

    seqs = []
    seq1_line = l.readline().strip()
    seq1_filename, seq1_start, seq1_end, seq1_strand, seq1_contig = seq1_line.strip().split(' ')
    seq1_filename = seq1_filename.replace("\"", "").rstrip('-')
    seqs.append((seq1_filename, seq1_start, seq1_end, seq1_strand, seq1_contig))

    seq2_line = l.readline().strip()
    seq2_filename, seq2_start, seq2_end, seq2_strand, seq2_contig = seq2_line.strip().split(' ')
    seq2_filename = seq2_filename.replace("\"", "").rstrip('-')
    seqs.append((seq2_filename,seq2_start, seq2_end, seq2_strand, seq2_contig))

    return seqs

def parse_h_stanza(l):
    ''' These are fasta headers, which give us coordinates for the window'''

    h_regex = re.compile("\"(\>(\w+\:\w+\-\w+)\"?(\s+\(reverse complement\))?\")")
    seq1_line = l.readline().strip()
    seq1_match = re.match(h_regex, seq1_line)
    if seq1_match.group(3) is None:
        sense = str("forward")
    else:
        sense = str("reverse")
    seq1_coords = seq1_match.group(2)

    seq2_line = l.readline().strip()
    seq2_match = re.match(h_regex, seq2_line)
    if seq2_match.group(3) is None:
        sense = str("forward")
    else:
        sense = str("reverse")
    seq2_coords = seq2_match.group(2)

    return ([(seq1_coords, sense), (seq2_coords, sense)])

def parse_a_stanza(l, linecount):
    '''These are the actual ungapped alignment blocks'''

    blockInfo = []
    alignment = []
    line = l.readline()
    block_len = None
    queryBlockStart = None
    subjBlockStart = None
    queryBlockEnd = None
    subjBlockEnd = None

    while line:
        fields = line.split()

        if fields[0] == "s":
            block_len = fields[1]
        elif fields[0] == "b":
            queryBlockStart = fields[1]
            subjBlockStart = fields[2]
        elif fields[0] == "e":
            queryBlockEnd = fields[1]
            subjBlockEnd = fields[2]
        elif fields[0] == "l":
            start1 = int(fields[1])
            start2 = int(fields[2])
            end1 = int(fields[3])
            end2 = int(fields[4])
            try:    pctId = int(fields[5])
            except: pctId = float(fields[5])

            if pctId >= 0.80:
                alignment.append((start1,end1,start2,end2,pctId))
        else:
            alnInfo = (block_len, queryBlockStart, queryBlockEnd, subjBlockStart, subjBlockEnd)
            blockInfo.append(alnInfo)
            blockInfo.append(alignment)

            return(blockInfo, linecount)

        linecount += 1
        line = l.readline()

def getOverlap(alnLoc, teLoc):
    """ If a TE overlaps this ungapped block then we want to return the coordinates for the overlap
        that correspond to the TE.
        + : what is returned
        - : Part of the alignment block that is not returned
        = : the TE
    """
    start = 0
    end = 0
    if max(0, min(alnLoc[1], teLoc[1]) - max(alnLoc[0], teLoc[0])) > 0:

        if alnLoc[0] <= teLoc[0] and alnLoc[1] >= teLoc[1]:
            #   ----+++++++----
            #       =======
            start = teLoc[0]
            end = teLoc[1]
        elif alnLoc[0] >= teLoc[0] and alnLoc[1] <= teLoc[1]:
            #     ++++++
            #  ===========
            start = alnLoc[0]
            end = alnLoc[1]
        elif alnLoc[1] >= teLoc[0] and alnLoc[1] <= teLoc[1]:
            #   ---++++++
            #      ======
            start = teLoc[0]
            end = alnLoc[1]
        elif alnLoc[0] >= teLoc[0] and alnLoc[0] <= teLoc[1]:
            #   +++++-----
            #   =====
            start = alnLoc[0]
            end = teLoc[1]
    else:
        start = None
        end = None
    return start,end

def merge_intervals(intervals):
    """ This will merge intervals that have any overlap into a single interval spanning both
    Help from here: https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    """
    # Sort the intervals in increasing order
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged

def getPercentCov(intervals, length):
    """ From all of the intervals overlapping the TE, calculate the percent coverage"""
    # Using merged_intervals calculate the total coverage for this element
    total_overlap = int(0)
    for interval in intervals:
        interval_start, interval_end = interval
        interval_length = interval_end - interval_start
        total_overlap += interval_length

    percentCov = round(total_overlap / length, 3)
    return percentCov

def getOverlaps(aStanza, queryElements, subjElements, queryCoords, subjCoords):
    """
    This subroutine is used to determine the amount of coverage
    of any element in the window.
    """
    # Regex that will be used to parse coordinates
    coord_regex = re.compile("(chr)?(\d+):(\d+)-(\d+)")

    # Now we can loop though all of the alignment blocks and see which overlap the TEs
    teOverlaps = defaultdict(lambda: defaultdict(list))
    for gapped_block in aStanza:
        blockInfo, alnInfo = gapped_block
        score, \
        queryBlockStart, \
        queryBlockEnd,  \
        subjBlockEndStart, \
        subjBlockEnd = blockInfo

        # These are the ungapped alignments
        for aln in alnInfo:
            gapQueryStart, \
            gapQueryEnd, \
            gapSubjStart, \
            gapSubjEnd, \
            perId = aln

            # Modify position in alignment to be the genomic position rather than alignment position
            modGapQueryStart = gapQueryStart + queryCoords[0]
            modGapQueryEnd = gapQueryEnd + queryCoords[0]
            modGapSubjStart = gapSubjStart + subjCoords[0]
            modGapSubjEnd = gapSubjEnd + subjCoords[0]

            # ASK: Do any of the TE's overlap this interval ?
            for element in queryElements:
                elementMatch = re.match(coord_regex, queryElements[element])
                elementStart = int(elementMatch.group(3))
                elementEnd = int(elementMatch.group(4))

                # Call overlap helper to actually calculate the overlap
                overlapStart, overlapEnd = getOverlap([modGapQueryStart, modGapQueryEnd], [elementStart, elementEnd])

                if overlapStart is not None and overlapEnd is not None:
                    teOverlaps[element]['query_coords'].append((overlapStart, overlapEnd))
                    modSubjStart = modGapSubjStart + (overlapStart - modGapQueryStart)
                    modSubjEnd = modGapSubjStart + (overlapEnd - modGapQueryStart)
                    teOverlaps[element]['subj_coords'].append((modSubjStart, modSubjEnd))

            '''Note: this should be cleaned up, is ugly'''
            # Repeat this for the subject elements
            for element in subjElements:
                elementMatch = re.match(coord_regex, subjElements[element])
                elementStart = int(elementMatch.group(3))
                elementEnd = int(elementMatch.group(4))
                # Call overlap helper to actually calculate the overlap
                overlapStart, overlapEnd = getOverlap([modGapSubjStart, modGapSubjEnd], [elementStart, elementEnd])
                if overlapStart is not None and overlapEnd is not None:
                    teOverlaps[element]['query_coords'].append((overlapStart, overlapEnd))
                    modQueryStart = modGapQueryStart + (overlapStart - modGapSubjStart)
                    modQueryEnd = modGapQueryStart + (overlapEnd - modGapSubjStart)
                    teOverlaps[element]['subj_coords'].append((modQueryStart, modQueryEnd))

    for te in teOverlaps:
        # If there was overlap, merge all of the intervals for this TE to prevent double-counting of any overlap
        mergedQueryIntervals = merge_intervals(teOverlaps[te]['query_coords'])
        mergedSubjIntervals = merge_intervals(teOverlaps[te]['subj_coords'])
        teOverlaps[te]['query_coords'] = mergedQueryIntervals
        teOverlaps[te]['subj_coords'] = mergedSubjIntervals

    return teOverlaps

def trimBlocks(aStanza):
    """ Trim ungapped blocks that are overlapping one another.
        Always use blocks from higher scoring stanza first and trim the blocks from the lower scoring stanza.
    For example:
      a {
        l 7941 8181 7961 8201 95
      }
      a {
        l 7955 12628 7990 12663 97
      }
      Gets trimmed to 7961 - 7990
      The subject sequence would also be trimmed by this same amount
    """

    alnLookup = defaultdict(list)
    trimmedAlns = []
    # First we create a lookup of all intervals
    for gapped_block in aStanza:
        blockInfo, alnInfo = gapped_block
        score = blockInfo[0]

        newAlns = []
        # Skip if this block was empty (no > 80% alignments)
        if not alnInfo:
            continue
        for aln in alnInfo:
            gapQueryStart, \
            gapQueryEnd, \
            gapSubjStart, \
            gapSubjEnd, \
            perId = aln
            # If this is the first ungapped alignment of the first gapped block
            if not alnLookup["query"] or not alnLookup["subj"]:
                # No need to check for overlap, we just store everything
                alnLookup["query"].append((gapQueryStart, gapQueryEnd))
                alnLookup["subj"].append((gapSubjStart, gapSubjEnd))
                newAlns.append((gapQueryStart, gapQueryEnd, gapSubjStart, gapSubjEnd, perId))
                queryBlockStart = int(gapQueryStart)
                subjBlockStart = int(gapSubjStart)
                queryBlockEnd = int(gapQueryEnd)
                subjBlockEnd = int(gapSubjEnd)
                continue

            # Does this interval overlap any previously stored interval
            for aln in alnLookup["query"]:
                overlapStart, overlapEnd = getOverlap([gapQueryStart, gapQueryEnd], aln)
                queryModLen = 0
                if overlapStart is not None or overlapEnd is not None:
                    newStart = overlapEnd
                    # If there was a complete overlap
                    queryModLen = int(newStart) - int(gapQueryStart)

            for aln in alnLookup["subj"]:
                overlapStart, overlapEnd = getOverlap([gapSubjStart, gapSubjEnd], aln)
                subjModLen = 0
                if overlapStart is not None or overlapEnd is not None:
                    newStart = overlapEnd
                    # If there was a complete overlap
                    subjModLen = int(newStart) - int(gapSubjStart)

            # Get the length to trimed based on most overlap between query/subj
            maxModLen = max(queryModLen, subjModLen)
            # Now adjust start coordinates by this length
            modQueryStart = int(gapQueryStart) + int(maxModLen)
            modSubjStart = int(gapSubjStart) + int(maxModLen)

            # If the full length was overlapping then remove it completely
            if int(modQueryStart) == int(gapQueryEnd) or int(modSubjStart) == int(gapSubjEnd):
                continue

            # Append this trimmed alignment to the dict
            alnLookup["query"].append((modQueryStart, gapQueryEnd))
            alnLookup["subj"].append((modSubjStart, gapSubjEnd))

            # If this is the first HSP of the ungapped block then store default start, end values
            if not newAlns:
                queryBlockStart = modQueryStart
                subjBlockStart = modSubjStart
                queryBlockEnd = gapQueryEnd
                subjBlockEnd = gapSubjEnd

            # Also append to list of new alignment blocks
            newAlns.append((modQueryStart, gapQueryEnd, modSubjStart, gapSubjEnd, perId))

            # Update the block endings
            queryBlockEnd = gapQueryEnd
            subjBlockEnd = gapSubjEnd

        # Rebuild the a-stanza structure
        newBlockInfo = (score, queryBlockStart, queryBlockEnd, subjBlockStart, subjBlockEnd)
        newAlnInfo = newAlns
        trimmedAlns.append([newBlockInfo, newAlnInfo])

    return(trimmedAlns)


def spit(sStanza, hStanza, aStanza, teLookup):
    window_name = []

    # Get some basic info on the flanks and their sense from the s-stanza
    queryFlank = sStanza[0][0]
    querySense = sStanza[0][3]
    subjFlank = sStanza[1][0]
    subjSense = sStanza[1][3]

    # Skipping reverse complement for now.
    if subjSense == "1":
        return

    # The window file stores windows ids like W22-W22_B73-B73, so we must recreate that string
    queryKey = queryFlank + "_" + subjFlank
    subjKey = subjFlank + "_" + queryFlank

    # We use this key to fetch all of the elements in this window
    queryElements = teLookup[queryKey]
    subjElements = teLookup[subjKey]

    # Regex that will be used to parse coordinates
    coord_regex = re.compile("(chr)?(\d+):(\d+)-(\d+)")

    # Get coordinate information for the query window
    queryCoord, queryStrand = hStanza[0]
    queryCoordMatch = re.match(coord_regex, queryCoord)
    queryStart = int(queryCoordMatch.group(3))
    queryEnd = int(queryCoordMatch.group(4))

    # Get coordinate info for the subject window
    subjCoord, subjStrand = hStanza[1]
    subjCoordMatch = re.match(coord_regex, subjCoord)
    subjStart = int(subjCoordMatch.group(3))
    subjEnd = int(subjCoordMatch.group(4))

    # We first, trim any ungapped alignment blocks that are overlapping with another block from a higher scoring block
    aStanzaTrim = trimBlocks(aStanza)

    # For each feature in this window, we compute it's coverage in the alignment using this subroutine.
    teOverlaps = getOverlaps(aStanzaTrim, queryElements, subjElements, [queryStart, queryEnd], [subjStart, subjEnd])

    #Now for every TE in this window, get the percent coverage and also determine which elements from the opposite genotype it hits to
    matches = defaultdict(lambda: defaultdict(dict))
    printStatements = []
    for queryTE in queryElements:
        # First we need to get the total length of each element to calculate percent coverage
        teCoords = re.match(coord_regex, queryElements[queryTE])
        teStart = int(teCoords.group(3))
        teEnd = int(teCoords.group(4))
        teLen = teEnd - teStart

        # Did any of the gapped alignments overlap this TE ?
        if queryTE in teOverlaps:
            percentCov = getPercentCov(teOverlaps[queryTE]['query_coords'], teLen)
        else:
            toPrint = "\t".join([str(queryTE), str(0), str("NA"), str("NA"), str(queryFlank), str("NA"), str(subjFlank), str("NA")])
            printStatements.append(toPrint)
            continue

        for subjTE in subjElements:
            teCoords = re.match(coord_regex, subjElements[subjTE])
            teStart = int(teCoords.group(3))
            teEnd = int(teCoords.group(4))
            teLen = teEnd - teStart

            for subjInterval in teOverlaps[queryTE]['subj_coords']:
                # Return the amount of overlap for this ungapped block and TE
                overlap = max(0, min(subjInterval[1], teEnd) - max(subjInterval[0], teStart))
                # If there was some overlap append the length of the overlap to the lookup
                if overlap > 0:
                    if subjTE not in matches[queryTE]:
                        matches[queryTE][subjTE] = overlap
                    else:
                        matches[queryTE][subjTE] += overlap
        # Now determine which elements mapped together the most and print out other relevant info
        bestHit = "NA"
        bestHitLen = int(0)
        match_list = []
        for match in matches[queryTE]:
            if match is None:
                match_list = "NA"
            match_list.append(match)
            if matches[queryTE][match] > bestHitLen:
                bestHit = match
                bestHitLen = matches[queryTE][match]

        query_coord = "-".join([ str(teOverlaps[queryTE]['query_coords'][0][0]), str(teOverlaps[queryTE]['query_coords'][-1][1]) ])
        subj_coord = "-".join([ str(teOverlaps[queryTE]['subj_coords'][0][0]), str(teOverlaps[queryTE]['subj_coords'][-1][1]) ])
        matches_to_print = ",".join(match_list)
        toPrint = "\t".join([str(queryTE), str(percentCov), str(bestHit), str(matches_to_print), str(queryFlank), str(query_coord), str(subjFlank), str(subj_coord)])
        printStatements.append(toPrint)

    # Repeat this process for the subject element NOTE: This should really be a subroutine of some kind, it's quite ineffiecient written this way.
    for subjTE in subjElements:
        # First we need to get the total length of each element to calculate percent coverage
        teCoords = re.match(coord_regex, subjElements[subjTE])
        teStart = int(teCoords.group(3))
        teEnd = int(teCoords.group(4))
        teLen = teEnd - teStart

        # Did any of the gapped alignments overlap this TE ?
        if subjTE in teOverlaps:
            percentCov = getPercentCov(teOverlaps[subjTE]['query_coords'], teLen)
        else:
            toPrint = "\t".join([str(subjTE), str(0), str("NA"), str("NA"), str(subjFlank), str("NA"), str(queryFlank), str("NA")])
            printStatements.append(toPrint)
            continue

        for queryTE in queryElements:
            teCoords = re.match(coord_regex, queryElements[queryTE])
            teStart = int(teCoords.group(3))
            teEnd = int(teCoords.group(4))
            teLen = teEnd - teStart

            for subjInterval in teOverlaps[subjTE]['subj_coords']:
                overlap = max(0, min(subjInterval[1], teEnd) - max(subjInterval[0], teStart))
                if overlap > 0:
                    if queryTE not in matches[subjTE]:
                        matches[subjTE][queryTE] = overlap
                    else:
                        matches[subjTE][queryTE] += overlap
        bestHit = "NA"
        bestHitLen = int(0)
        match_list = []
        for match in matches[subjTE]:
            if match is None:
                match_list = "NA"
            match_list.append(match)
            if matches[subjTE][match] > bestHitLen:
                bestHit = match
                bestHitLen = matches[subjTE][match]

        query_coord = "-".join([ str(teOverlaps[subjTE]['query_coords'][0][0]), str(teOverlaps[subjTE]['query_coords'][-1][1]) ])
        subj_coord = "-".join([ str(teOverlaps[subjTE]['subj_coords'][0][0]), str(teOverlaps[subjTE]['subj_coords'][-1][1]) ])
        matches_to_print = ",".join(match_list)
        toPrint = "\t".join([str(subjTE), str(percentCov), str(bestHit), str(matches_to_print), str(subjFlank), str(query_coord), str(queryFlank), str(subj_coord)])
        printStatements.append(toPrint)

    return printStatements

def storeFeatures(queryFile, subjFile):
    """
    Loop through the window files and store information on every TE using the flanking genes as a key in the dictionary.
    """

    teLookup = defaultdict(lambda: defaultdict(dict))
    for elementFile in [queryFile, subjFile]:
        with open(elementFile) as f:
            for line in f:
                element, element_loc, queryFlank, queryFlankLoc, subjFlank, subjFlankLoc = line.strip().split("\t")
                key = "_".join([queryFlank,subjFlank])
                teLookup[key][element] = element_loc

    return teLookup

def parselav(queryTE, subjTE, lavfile):
    """Read through the LAV File and get info from the different stanzas to pass to spit subroute that
       "spits" out info on how well each element aligned.
    """
    s_info = None
    h_info = None
    block = None
    alnBlocks = []
    teHash = storeFeatures(queryTE, subjTE)

    with open(lavfile) as l:
        line = l.readline()
        count = int(0)
        while line:
            count += 1
            line = l.readline()

            if line.startswith("#:lav") or line.startswith("#:eof"):
                if s_info is None:
                    continue
                # I am not using the reverse complement alignmnets for now.
                if h_info[0][1] == "reverse" or h_info[1][1] == "reverse":
                    continue
                # This subroutine will parse all of the stored a stanzas.
                lines = spit(s_info, h_info, alnBlocks, teHash)
                for line in lines:
                    print(line)
                continue
            elif line.startswith("d {"):
                parse_d_stanza(l)
            elif line.startswith("s {"):
                s_info = parse_s_stanza(l)
            elif line.startswith("h {"):
                h_info = parse_h_stanza(l)
            elif line.startswith("a {"):
                block, count = parse_a_stanza(l, count)
                alnBlocks.append(block)
            else:
                continue


def main():
    '''Main subroutine that makes calls to all relevant functions'''

    if len(sys.argv) == 4:
        for filename in glob.glob(os.path.join(sys.argv[1], '*.lav')):
            print(filename)
            parselav(sys.argv[2], sys.argv[3], filename)

if __name__ == '__main__':
    main()
