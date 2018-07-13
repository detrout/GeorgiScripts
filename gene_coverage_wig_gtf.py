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
    parser.add_argument('-field1', dest='sourcetype', help='sourcetype name')
    parser.add_argument('-normalize', action='store_true', default=False,
                        help='normalize scores')
    parser.add_argument('-maxGeneLength', default=None, type=int,
                        help='maximum gene length to consider')
    parser.add_argument('-singlemodelgenes', default=False, action='store_true',
                        help='only consider genes with a single model')
    parser.add_argument('-genetype', help='limit to specified gene type')
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

    if args.sourcetype is not None:
        logger.info('will only consider %s genes', args.sourceType)

    if args.maxGeneLength is not None:
        logger.info('will only consider genes longer than %s and shorter than %s',
                    args.minGeneLength,
                    args.maxGeneLength)

    if args.normalize:
        logger.info('will normalize scores')

    if args.singlemodelgenes:
        logger.info('will only use genes with one isoform')

    if args.genetype:
        logger.info('will only consider genes of type %s', args.genetype)

    GeneDict = get_gene_dict(args.gtf, args.sourcetype, args.genetype)

    CoverageDict = build_coverage_dict(GeneDict, args.singlemodelgenes)
    score_wiggle(args.wig, CoverageDict)

    if args.printlist:
        outfile = open(args.outputfilename + '.geneList', 'wt')
    else:
        outfile = None

    outputArray, geneNumber = compute_coverage_array(
            GeneDict, CoverageDict,
            args.minGeneLength, args.maxGeneLength,
            outfile)
    if args.printlist:
        outfile.close()

    outfile=open(args.outputfilename, 'wt')

    for i in range(100):
        if args.normalize:
            outfile.write(str(i) + '\t' + str(outputArray[i]/geneNumber)+'\n')
        else:
            outfile.write(str(i) + '\t' + str(outputArray[i])+'\n')

    outfile.close()

def get_gene_dict(filename, source, gene_type_filter=None):
    listoflines = open(filename, 'rt')
    GeneDict={}
    for line in listoflines:
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
        if geneID in GeneDict:
            pass
        else:
            GeneDict[geneID]={}
        if transcriptID in GeneDict[geneID]:
            if chromosome != GeneDict[geneID][transcriptID][0][0]:
                continue
        else:
            GeneDict[geneID][transcriptID]=[]
        GeneDict[geneID][transcriptID].append((chromosome,left,right,strand))

    logger.info('finished inputting annotation %s', len(geneDict.keys()))
    return GeneDict

def get_gff_attribute_value_by_key(field, name):
    start = field.index(name)
    start += len(name)
    start = field.index('"', start)
    start += 1
    end = field.index('"', start)
    return field[start:end]

def build_coverage_dict(GeneDict, singlemodelgenes):
    CoverageDict = {}
    genesToRemove = set()
    for geneID in GeneDict:
        if singlemodelgenes and len(GeneDict[geneID]) > 1:
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

def score_wiggle(wigglename, CoverageDict):
    listoflines = open(wigglename, 'rt')
    for line in listoflines:
        if line.startswith('track'):
            continue
        if line.startswith('#'):
            continue
        fields=line.replace(' ','\t').strip().split('\t')
        chromosome=fields[0]
        left=int(fields[1])
        right=int(fields[2])
        score=float(fields[3])
        if chromosome in CoverageDict:
            pass
        else:
            continue
        for j in range(left,right):
            if j in CoverageDict[chromosome]:
                CoverageDict[chromosome][j]=score

    logger.info('finished inputting wiggle')
    logger.info('genes passed type filters %s', len(geneDict))

def compute_coverage_array(GeneDict, CoverageDict, minGeneLength, maxGeneLength, geneStream=None):
    outputArray = numpy.zeros(shape=100)

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
                counts.append(CoverageDict[chromosome][NucleotideList[i]])
            final_vector.append(numpy.mean(counts))
            k=b
        i=0
        for v in final_vector:
            outputArray[i]+=v
            i+=1
        geneNumber+=1
        if geneStream is not None:
            geneStream.write(geneID)
            geneStream.write('\n')

    logger.info('%s genes considered', geneNumber)
    return outputArray, geneNumber

if __name__ == '__main__':
    main()
