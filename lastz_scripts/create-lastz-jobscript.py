#!/usr/bin/env python3

import sys
from collections import defaultdict

def jobs(window_file):
    bashShebang = "#!/usr/bin/bash"
    modules = ["module load lastz/1.04.00", "module load samtools"]
    joboptions = ["#PBS -l walltime=20:00:00,nodes=1:ppn=1,mem=8gb", "#PBS -V", "#PBS -N run-lastz"]
    programoptions = [ "--format=lav", "--gfextend", "--chain", "--ambiguous=n", "--gap=450,30", "--gapped", "--markend", "--hspthresh=3500", "--gappedthresh=3500", "--filter=identity:80", "--allocate:traceback=524288000"]
    workdir = ["/scratch.global/broha006/projects/w22/lastz_outputs/w22_ph207"]
    queryGenome = "/home/maize/shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta"
    subjGenome = "/home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa"

    windows = defaultdict(lambda: defaultdict(dict))
    with open(window_file) as w:
        for line in w:
            element, elementLoc, query, queryLoc, subj, subjLoc = line.strip().split('\t')
            if "ctg" in queryLoc: 
                 continue
            if "ctg" in subjLoc:
                 continue 
            if "scaffold" in subjLoc:
                 continue
            if queryLoc == "NA" or subjLoc == "NA":
                 continue
            windows[queryLoc][subjLoc] = "_".join([query,subj])

    print(bashShebang)
    for joboption in joboptions:
        print(joboption)
    for module in modules:
        print(module)
    print("cd " + str(workdir[0]))

    for queryCoord in windows:
        for subjCoord in windows[queryCoord]:
            queryOut, subjOut = windows[queryCoord][subjCoord].split('_')
            querySam = " ".join(["samtools faidx", queryGenome, "chr" + queryCoord, ">", queryOut])
            subjSam = " ".join(["samtools faidx", subjGenome, "chr0" + subjCoord, ">", subjOut])
            opts = " ".join(programoptions)
            runAln = " ".join(["lastz", subjOut, queryOut, opts + " --rdotplot=" + str(windows[queryCoord][subjCoord]) + ".txt", ">", str(windows[queryCoord][subjCoord]) + ".lav"])
            cleanup = "".join(["rm ", queryOut + ";", "rm ", subjOut])
            command = ";".join([querySam, subjSam, runAln, cleanup])
            print(command)

def main():
    '''Main subroutine that makes calls to all relevant functions'''

    if len(sys.argv) != 2:
        print("Need 1 arguments")
    else:
        jobs(sys.argv[1])

if __name__ == '__main__':
	main()
