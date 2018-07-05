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

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-field1', dest='sourcetype', help='sourcetype name')
    parser.add_argument('-normalize', action='store_true', default=False,
                        help='normalize scores')
    parser.add_argument('-maxGeneLength', default=None, type=int,
                        help='maximum gene length to consider')
    parser.add_argument('-singlemodelgenes', default=False, action='store_true',
                        help='only consider genes with a single model')
    parser.add_argument('-printlist', default=False, action='store_true',
                        help='write gene ids considered to <outfile>.geneList')
    parser.add_argument('gtf', help='GTF file name')
    parser.add_argument('wig', help='wiggle file to score')
    parser.add_argument('minGeneLength', type=int, help='minimum gene length to consider')
    parser.add_argument('outputfilename', help='filename to write coverage data to')
    args = parser.parse_args(cmdline)

    if args.sourcetype is not None:
        print('will only consider', args.sourcetype, 'genes')

    if args.maxGeneLength is not None:
        print('will only consider genes longer than', args.minGeneLength, 'and shorter than', args.maxGeneLength)

    if args.normalize:
        print('will normalize scores')

    if args.singlemodelgenes:
        print('will only use genes with one isoform')

    GeneDict = get_gene_dict(args.gtf, args.sourcetype)
    print('finished inputting annotation')

    CoverageDict = build_coverage_dict(GeneDict, args.singlemodelgenes)
    score_wiggle(args.wig, CoverageDict)
    print('finished inputting wiggle')

    print(len(GeneDict.keys()))

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

def get_gene_dict(filename, source):
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
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if GeneDict.has_key(geneID):
            pass
        else:
            GeneDict[geneID]={}
        if GeneDict[geneID].has_key(transcriptID):
            if chromosome != GeneDict[geneID][transcriptID][0][0]:
                continue
        else:
            GeneDict[geneID][transcriptID]=[]
        GeneDict[geneID][transcriptID].append((chromosome,left,right,strand))

    return GeneDict

def build_coverage_dict(GeneDict, singlemodelgenes):
    CoverageDict = {}
    i=0
    for geneID in GeneDict.keys():
        if singlemodelgenes and len(GeneDict[geneID].keys()) > 1:
            del GeneDict[geneID]
            continue
        i+=1
        for transcriptID in GeneDict[geneID].keys():
            for (chromosome,left,right,strand) in GeneDict[geneID][transcriptID]:
                if CoverageDict.has_key(chromosome):
                    pass
                else:
                    CoverageDict[chromosome]={}
                for j in range(left,right):
                    CoverageDict[chromosome][j]=0
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
        if CoverageDict.has_key(chromosome):
            pass
        else:
            continue
        for j in range(left,right):
            if CoverageDict[chromosome].has_key(j):
                CoverageDict[chromosome][j]=score


def compute_coverage_array(GeneDict, CoverageDict, minGeneLength, maxGeneLength, geneStream=None):
    outputArray = numpy.zeros(shape=100)

    geneNumber=0.0
    for geneID in GeneDict.keys():
        NucleotideList=[]
        for transcriptID in GeneDict[geneID].keys():
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

    print(geneNumber, 'genes considered')
    return outputArray, geneNumber

if __name__ == '__main__':
    main()
