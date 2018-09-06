# Initial Code By:
# Georgi Marinov 05/19/2012
#
# Version 1.0.1
from __future__ import print_function

import sys
import pysam
import argparse


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf", help="gtf file")
    parser.add_argument("SAM", help="SAM file")
    parser.add_argument("chrominfo", help="chrominfo file")
    parser.add_argument("outputfilename", help="output file name")
    parser.add_argument("--nomulti", help="discard multireads",
                        action="store_true", dest="noMulti")
    parser.add_argument("--nounique", action="store_true", dest="noUnique",
                        help="discard unique reads")
    parser.add_argument("--noNH", action="store_true", dest="noNH",
                        help="use if multi reads do not have NH tags",)
    parser.add_argument("--verbose", action="store_true",
                        help="verbose output")

    return parser


def main(cmdline=None):
    # maybe bedtools can just replace all of this
    # also see samtools idx stats
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.verbose:
        if args.noMulti:
            print('will discard multi reads')
        if args.noUnique:
            print('will discard unique reads')

    with open(args.chrominfo) as chromData:
        chromInfoList = getChromInfo(chromData)

    samfile = pysam.Samfile(args.SAM, "rb")
    if args.noNH:
        ReadMultiplicityDict = getReadMultiplicity(
            chromInfoList, samfile, verbose=args.verbose)
    else:
        ReadMultiplicityDict = {}

    PosCountsDict, TotalReads = pileUpReads(
        chromInfoList, samfile, ReadMultiplicityDict, noNH=False,
        noMulti=False, verbose=args.verbos)
    if args.verbose:
        print('....................................')

    with open(args.gtf) as gtfData:
        geneDict = readGtf(gtfData)
    if args.verbose:
        print('finished inputting annotation')

    ExonicReads, IntronicReads = getCounts(PosCountsDict, geneDict)
    with open(args.outfilename, 'w') as outfile:
        writeOutput(outfile, TotalReads, ExonicReads, IntronicReads)


def getChromInfo(chrominfo):
    chromInfoList = []
    for line in chrominfo:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = 0
        end = int(fields[1])
        chromInfoList.append((chrom, start, end))

    return chromInfoList


def getReadMultiplicity(chromInfoList, samfile, verbose=False):
    readMultiplicityDict = {}
    i = 0
    for (chrom, start, end) in chromInfoList:
        # this just checks to see if there is an entry by looking
        # instead of error handling
        try:
            jj = 0
            for alignedread in samfile.fetch(chrom, start, end):
                jj += 1
                if jj == 1:
                    break
        except ValueError:
            if verbose:
                print(chrom, start, end, 'not found in BAM file, skipping')
            continue

        # here is where the work begins
        for alignedread in samfile.fetch(chrom, start, end):
            if verbose:
                i += 1
                if i % 5000000 == 0:
                    print(
                        str(i/1000000) +
                        'M alignments processed in multiplicity examination',
                        chrom, start, alignedread.pos, end)

            fields = str(alignedread).split('\t')
            ID = fields[0]
            if alignedread.is_read1:
                ID = ID + '/1'

            if alignedread.is_read2:
                ID = ID + '/2'

            try:
                readMultiplicityDict[ID] += 1
            except KeyError:
                readMultiplicityDict[ID] = 1

    return readMultiplicityDict


def pileUpReads(chromInfoList, samfile, readMultiplicityDict,
                noNH=False, noMulti=False, verbose=False):
    posCountsDict = {}
    i = 0
    totalReads = 0.0
    for (chrom, start, end) in chromInfoList:
        # same checking code here
        try:
            jj = 0
            for alignedread in samfile.fetch(chrom, start, end):
                jj += 1
                if jj == 1:
                    break
        except ValueError:
            if verbose:
                print('{} {} {} not found in BAM file, skipping'.format(
                    chrom, start, end))
            continue

        if chrom in posCountsDict:
            pass
        else:
            posCountsDict[chrom] = {}

        for alignedread in samfile.fetch(chrom, start, end):
            if verbose:
                i += 1
                if i % 5000000 == 0:
                    print('{} M alignments processed'.format(i/1000000))

            if noNH:
                fields = str(alignedread).split('\t')
                ID = fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'

                if alignedread.is_read2:
                    ID = ID + '/2'

                scaleby = readMultiplicityDict[ID]
            else:
                try:
                    scaleby = alignedread.get_tag('NH')
                except KeyError:
                    print('multireads not specified with the NH tag, exiting')
                    sys.exit(1)

            if noMulti and scaleby > 1:
                continue

            weight = 1.0/scaleby
            totalReads += weight
            pos = alignedread.pos
            try:
                posCountsDict[chrom][pos] += weight
            except KeyError:
                posCountsDict[chrom][pos] = weight

    return posCountsDict, totalReads


def getCounts(posCountsDict, geneDict, verbose=False):
    exonicReads = 0
    intronicReads = 0
    keys = sorted(geneDict.keys())
    for chrom in keys:
        if verbose:
            print(chrom)

        for geneID in geneDict[chrom].keys():
            coordinates = []
            for (start, stop) in geneDict[chrom][geneID]:
                coordinates.append(start)
                coordinates.append(stop)
                for i in range(start, stop):
                    try:
                        exonicReads += posCountsDict[chrom][i]
                        # delete entry so that it is not counted multiple times
                        # either because of multiple models or the next
                        # intron pass
                        del posCountsDict[chrom][i]
                    except KeyError:
                        pass

            # all the exonic reads have been counted and removed
            # so the remainder in the rage is intronic
            for i in range(min(coordinates), max(coordinates)):
                try:
                    intronicReads += posCountsDict[chrom][i]
                    del posCountsDict[chrom][i]
                except KeyError:
                    pass

    return exonicReads, intronicReads


def readGtf(gtf):
    geneDict = {}
    for line in gtf:
        if line[0] == '#':
            continue

        fields = line.strip().split('\t')
        if fields[2] != 'exon':
            continue

        chrom = fields[0]
        if chrom in geneDict:
            pass
        else:
            geneDict[chrom] = {}

        start = int(fields[3])
        stop = int(fields[4])
        geneID = fields[8].split('gene_id "')[1].split('";')[0]

        try:
            geneDict[chrom][geneID].add((start, stop))
        except KeyError:
            geneDict[chrom][geneID] = set([(start, stop)])

    return geneDict


def writeOutput(outfile, totalReads, exonicReads, intronicReads):
    outfile.write('#Class\tFraction\n')
    outfile.write('Exonic:{}\t\n'.format(exonicReads/totalReads))
    outfile.write('Intronic:{}\t\n'.format(intronicReads/totalReads))
    intergenicReads = totalReads - exonicReads - intronicReads
    outfile.write('Intergenic:{}\t\n'.format(intergenicReads/totalReads))


class pseudoSAM():
    def __init__(self, readData):
        self.reads = []

    def fetch(self, chrom, start, stop):
        return self.reads


class pseudoFile():
    def __init__(self, readLineData):
        self.readLineEntries = readLineData

    def write(self, dataEntry):
        self.readLineEntries.append(dataEntry)

    def read(self):
        return self.readLineEntries.pop()


if __name__ == '__main__':
    main()
