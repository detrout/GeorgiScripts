##################################
#                                #
# Last modified 04/05/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################
from __future__ import print_function

from argparse import ArgumentParser
import sys
import numpy
import logging

logger = logging.getLogger('gene_coverage_wig_gtf')


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-field1', dest='sourceType', help='sourcetype name')
    parser.add_argument('-normalize', action='store_true', default=False,
                        help='normalize scores')
    parser.add_argument('-maxGeneLength', default=None, type=int,
                        help='maximum gene length to consider')
    parser.add_argument('-singlemodelgenes', default=False, action='store_true',
                        dest='singleModelGenes',
                        help='only consider genes with a single model')
    parser.add_argument('-genetype', dest='geneType', help='limit to specified gene type')
    parser.add_argument('-printlist', default=False, action='store_true',
                        help='write gene ids considered to <outfile>.geneList')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='only report errors')
    parser.add_argument('gtf', help='GTF file name')
    parser.add_argument('wig', help='wiggle file to score')
    parser.add_argument('minGeneLength', type=int, help='minimum gene length to consider')
    parser.add_argument('outputfilename', help='filename to write coverage data to')
    args = parser.parse_args(cmdline)

    if args.quiet:
        logging.basicConfig(level=logging.WARN)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.sourceType is not None:
        logger.info('will only consider %s genes', args.sourceType)

    if args.maxGeneLength is not None:
        logger.info('will only consider genes longer than %s and shorter than %s',
                    args.minGeneLength,
                    args.maxGeneLength)

    if args.normalize:
        logger.info('will normalize scores')

    if args.singleModelGenes:
        logger.info('will only use genes with one isoform')

    if args.geneType:
        logger.info('will only consider genes of type %s', args.genetype)

    with open(args.gtf) as gtfStream:
        geneDict = readAnnotation(
            gtfStream,
            args.sourceType,
            args.geneType)

    with open(args.wig) as wigStream:
        coverageDict = readWiggle(wigStream, geneDict, args.singleModelGenes)

    if args.printlist:
        geneListFilename = args.outputfilename + '.geneList'
    else:
        geneListFilename = None

    outputArray = createCoverageArray(
        geneDict, coverageDict,
        args.minGeneLength, args.maxGeneLength,
        geneListFilename,
        args.normalize
    )

    with open(args.outputfilename, 'wt') as outfile:
        for i in range(100):
            outfile.write(str(i) + '\t' + str(outputArray[i])+'\n')

def readAnnotation(stream, source, doSingleModel, gene_type_filter=None):
    geneDict={}
    entriesToDelete=set()
    for line in stream:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if source is not None and fields[1] != source:
            continue
        chromosome=fields[0]
        strand=fields[6]
        left=int(fields[3])
        right=int(fields[4])
        if gene_type_filter is not None:
            try:
                geneType = get_gff_attribute_value_by_key(fields[8], 'gene_type')
            except ValueError:
                continue
            if geneType != gene_type_filter:
                continue
        geneID = get_gff_attribute_value_by_key(fields[8], 'gene_id')
        transcriptID = get_gff_attribute_value_by_key(fields[8], 'transcript_id')
        if geneID in geneDict:
            pass
        else:
            geneDict[geneID]={}
        if transcriptID in geneDict[geneID]:
            if chromosome != geneDict[geneID][transcriptID][0][0]:
                continue
        else:
            geneDict[geneID][transcriptID]=[]
        geneDict[geneID][transcriptID].append((chromosome,left,right,strand))

    logger.info('finished inputting annotation %s', len(geneDict.keys()))
    return geneDict

def get_gff_attribute_value_by_key(field, name):
    start = field.index(name)
    start += len(name)
    start = field.index('"', start)
    start += 1
    end = field.index('"', start)
    return field[start:end]

def initializeCoverageDict(GeneDict, singleModelGenes):
    CoverageDict = {}
    genesToRemove = set()
    for geneID in GeneDict:
        if singleModelGenes and len(GeneDict[geneID]) > 1:
            genesToRemove.add(geneID)
            continue
        for transcriptID in GeneDict[geneID]:
            for (chromosome,left,right,strand) in GeneDict[geneID][transcriptID]:
                if chromosome in CoverageDict:
                    pass
                else:
                    CoverageDict[chromosome]={}
                for j in range(left,right):
                    CoverageDict[chromosome][j]=0
    for geneID in genesToRemove:
        del GeneDict[geneID]
    return CoverageDict

def readWiggle(wiggle, geneDict, singleModelGenes):
    coverageDict = initializeCoverageDict(geneDict, singleModelGenes)
    for line in wiggle:
        if line.startswith('track'):
            continue
        if line.startswith('#'):
            continue
        fields=line.replace(' ','\t').strip().split('\t')
        chromosome=fields[0]
        left=int(fields[1])
        right=int(fields[2])
        score=float(fields[3])
        if chromosome in coverageDict:
            pass
        else:
            continue
        for j in range(left,right):
            if j in coverageDict[chromosome]:
                coverageDict[chromosome][j]=score

    logger.info('finished inputting wiggle')
    logger.info('genes passed type filters %s', len(geneDict))
    return coverageDict

def createCoverageArray(GeneDict, coverageDict,
                        minGeneLength, maxGeneLength=None,
                        geneListFilename=None,
                        doNormalize=False):
    outputArray = numpy.zeros(shape=100)

    geneListStream = open(geneListFilename, 'wt') if geneListFilename else None

    geneNumber=0.0
    for geneID in GeneDict:
        NucleotideList=[]
        for transcriptID in GeneDict[geneID]:
            for (chromosome,left,right,strand) in GeneDict[geneID][transcriptID]:
                for i in range(left,right):
                    NucleotideList.append(i)
        NucleotideList=sorted(set(NucleotideList))
        if strand=='-' or strand=='R':
            NucleotideList.reverse()
        geneLength=len(NucleotideList)
        if geneLength < minGeneLength:
            continue
        if maxGeneLength is not None and geneLength > maxGeneLength:
            continue
        stepsize = geneLength/100.0
        k=0
        final_vector=[]
        while k < geneLength-stepsize:
            b=k+stepsize
            counts=[]
            for i in range(int(k),int(b)):
                counts.append(coverageDict[chromosome][NucleotideList[i]])
            final_vector.append(numpy.mean(counts))
            k=b
        i=0
        for v in final_vector:
            outputArray[i]+=v
            i+=1
        geneNumber+=1
        if geneListStream:
            geneListStream.write(geneID)
            geneListStream.write('\n')

    logger.info('%s genes considered', geneNumber)

    if doNormalize:
        outputArray /= geneNumber

    if geneListStream:
        geneListStream.close()

    return outputArray


if __name__ == '__main__':
    main()
