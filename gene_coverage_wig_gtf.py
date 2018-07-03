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
from sets import Set

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-field1', dest='biotype', help='biotype name')
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

    if args.biotype is not None:
        print('will only consider', args.biotype, 'genes')

    if args.maxGeneLength is not None:
        print('will only consider genes longer than', args.minGeneLength, 'and shorter than', args.maxGeneLength)

    if args.normalize:
        print('will normalize scores')

    if args.singlemodelgenes:
        print('will only use genes with one isoform')

    listoflines = open(args.gtf, 'rt')
    GeneDict={}
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if args.biotype is not None:
            if fields[1] != args.biotype:
                continue
        chr=fields[0]
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
            if chr != GeneDict[geneID][transcriptID][0][0]:
                continue
        else:
            GeneDict[geneID][transcriptID]=[]
        GeneDict[geneID][transcriptID].append((chr,left,right,strand))

    print('finished inputting annotation')

    CoverageDict={}

    i=0
    for geneID in GeneDict.keys():
        if args.singlemodelgenes and len(GeneDict[geneID].keys()) > 1:
            del GeneDict[geneID]
            continue
        i+=1
        for transcriptID in GeneDict[geneID].keys():
            for (chr,left,right,strand) in GeneDict[geneID][transcriptID]:
                if CoverageDict.has_key(chr):
                    pass
                else:
                    CoverageDict[chr]={}
                for j in range(left,right):
                    CoverageDict[chr][j]=0

    listoflines = open(args.wig, 'rt')
    for line in listoflines:
        if line.startswith('track'):
            continue
        if line.startswith('#'):
            continue
        fields=line.replace(' ','\t').strip().split('\t')
        chr=fields[0]
        left=int(fields[1])
        right=int(fields[2])
        score=float(fields[3])
        if CoverageDict.has_key(chr):
            pass
        else:
            continue
        for j in range(left,right):
            if CoverageDict[chr].has_key(j):
                CoverageDict[chr][j]=score

    print('finished inputting wiggle')

    output_Array={}
    for i in range(100):
        output_Array[i]=0

    print(len(GeneDict.keys()))

    if args.printlist:
        outfile=open(args.outputfilename + '.geneList','w')

    geneNumber=0.0
    for geneID in GeneDict.keys():
        NucleotideList=[]
        for transcriptID in GeneDict[geneID].keys():
            for (chr,left,right,strand) in GeneDict[geneID][transcriptID]:
                for i in range(left,right):
                    NucleotideList.append(i)
        NucleotideList=list(Set(NucleotideList))
        NucleotideList.sort()
        if strand=='-' or strand=='R':
            NucleotideList.reverse()
        geneLength=len(NucleotideList)
        if geneLength < args.minGeneLength:
            continue
        if args.maxGeneLength is not None and geneLength > args.maxGeneLength:
            continue
        stepsize = geneLength/100.0
        k=0
        final_vector=[]
        while k < geneLength-stepsize:
            b=k+stepsize
            counts=[]
            for i in range(int(k),int(b)):
                counts.append(CoverageDict[chr][NucleotideList[i]])
            final_vector.append(numpy.mean(counts))
            k=b
        i=0
        for v in final_vector:
            output_Array[i]+=v
            i+=1
        geneNumber+=1
        if args.printlist:
            outfile.write(geneID + '\n')

    if args.printlist:
        outfile.close()

    print(geneNumber, 'genes considered')

    outfile=open(args.outputfilename, 'wt')

    for i in range(100):
        if args.normalize:
            outfile.write(str(i) + '\t' + str(output_Array[i]/geneNumber)+'\n')
        else:
            outfile.write(str(i) + '\t' + str(output_Array[i])+'\n')

    outfile.close()

if __name__ == '__main__':
    main()
