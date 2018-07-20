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
import unittest

logger = logging.getLogger('gene_coverage_wig_gtf')


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.quiet:
        logging.basicConfig(level=logging.WARN)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.source_type is not None:
        logger.info('will only consider genes from source %s', args.source_type)

    if args.max_gene_length is not None:
        logger.info('will only consider genes longer than %s and shorter than %s',
                    args.min_gene_length,
                    args.max_gene_length)
    else:
        logger.info('will only consder genes longer than %s', args.min_gene_length)

    if args.normalize:
        logger.info('will normalize scores')

    if args.all_gene_models:
        logger.info('will only all genes')
    else:
        logger.info('will only use genes with one isoform')

    if args.gene_type:
        logger.info('will only consider genes of type %s', args.gene_type)

    with open(args.gtf) as gtfStream:
        geneDict = readAnnotation(
            gtfStream,
            args.source_type,
            args.gene_type)

    with open(args.wig) as wigStream:
        coverageDict = readWiggle(wigStream, geneDict, args.all_gene_models)

    if args.print_list:
        geneListFilename = args.outputfilename + '.geneList'
    else:
        geneListFilename = None

    outputArray = createCoverageArray(
        geneDict, coverageDict,
        args.min_gene_length, args.max_gene_length,
        geneListFilename,
        args.normalize
    )

    with open(args.output, 'wt') as outfile:
        for i in range(100):
            outfile.write(str(i) + '\t' + str(outputArray[i])+'\n')

def make_parser():
    parser = ArgumentParser()
    parser.add_argument('--gtf', required=True, help='GTF file name')
    parser.add_argument('wig', help='wiggle file to score')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    parser.add_argument('--source-type', help='source type name')
    parser.add_argument('--max-gene-length', default=None, type=int,
                        help='maximum gene length to consider')
    parser.add_argument('--min-gene-length', type=int, default=1000,
                        help='minimum gene length to consider')
    parser.add_argument('--all-gene-models', default=False, action='store_true',
                        help='use all gene models, not just the ones with a single model')
    parser.add_argument('--gene-type', help='limit to specified gene type')
    parser.add_argument('--print-list', default=False, action='store_true',
                        help='write gene ids considered to <outfile>.geneList')
    parser.add_argument('--normalize', action='store_true', default=False,
                        help='normalize scores')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='only report errors')
    parser.add_argument('--version', action='store_true', help='report version number')
    return parser

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
                geneType = getGFFAttributeValueByKey(fields[8], 'gene_type')
            except ValueError:
                continue
            if geneType != gene_type_filter:
                continue
        geneID = getGFFAttributeValueByKey(fields[8], 'gene_id')
        transcriptID = getGFFAttributeValueByKey(fields[8], 'transcript_id')
        if transcriptID in geneDict.setdefault(geneID, {}):
            if chromosome != geneDict[geneID][transcriptID][0][0]:
                continue
        else:
            geneDict[geneID][transcriptID]=[]
        geneDict[geneID][transcriptID].append((chromosome,left,right,strand))

    logger.info('finished inputting annotation %s', len(geneDict.keys()))
    return geneDict

def getGFFAttributeValueByKey(field, name):
    start = field.index(name)
    start += len(name)
    start = field.index('"', start)
    start += 1
    end = field.index('"', start)
    return field[start:end]

def initializeCoverageDict(GeneDict, all_gene_models):
    CoverageDict = {}
    genesToRemove = set()
    for geneID in GeneDict:
        if all_gene_models == False and len(GeneDict[geneID]) > 1:
            genesToRemove.add(geneID)
            continue
        for transcriptID in GeneDict[geneID]:
            for (chromosome,left,right,strand) in GeneDict[geneID][transcriptID]:
                for j in range(left,right):
                    CoverageDict.setdefault(chromosome, {})[j]=0
    for geneID in genesToRemove:
        del GeneDict[geneID]
    return CoverageDict

def readWiggle(wiggle, geneDict, all_gene_models):
    coverageDict = initializeCoverageDict(geneDict, all_gene_models)
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
        if chromosome not in coverageDict:
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
        NucleotideList, chromosome = buildNucleotideList(GeneDict[geneID])
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
            geneListStream.write('\t')
            geneListStream.write('\t'.join([str(x) for x in final_vector]))
            geneListStream.write('\n')

    logger.info('%s genes considered', geneNumber)

    if doNormalize:
        outputArray /= geneNumber

    if geneListStream:
        geneListStream.close()

    return outputArray


def buildNucleotideList(geneModel):
    NucleotideList=[]
    for transcriptID in geneModel:
        for (chromosome,left,right,strand) in geneModel[transcriptID]:
            for i in range(left,right):
                NucleotideList.append(i)
    NucleotideList=sorted(set(NucleotideList))
    if strand=='-' or strand=='R':
        NucleotideList.reverse()
    return NucleotideList, chromosome


class GeneCoverageWigGtfTest(unittest.TestCase):
    def testInitializeCoverageDictEmpty(self):
        coverageDict = {}
        geneDict = {}
        all_gene_models = True
        self.assertEqual(coverageDict, initializeCoverageDict(geneDict, all_gene_models))

    def testInitializeCoverageDictSimple(self):
        coverageDict = {"chr1": {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")]}}
        all_gene_models = True
        self.assertEqual(coverageDict, initializeCoverageDict(geneDict, all_gene_models))

    def testInitializeCoverageDictMulti(self):
        coverageDict = {"chr1": {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0, 15: 0, 16: 0, 17: 0, 18: 0,
                                 19: 0},
                        "chr2": {20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}}
        all_gene_models = True
        self.assertEqual(coverageDict, initializeCoverageDict(geneDict, all_gene_models))

    def testReadWiggleEmpty(self):
        coverageDict = {}
        wig = []
        geneDict = {}
        all_gene_models = True
        self.assertEqual(coverageDict, readWiggle(wig, geneDict, all_gene_models))

    def testReadWiggleSingle(self):
        coverageDict = {"chr1": {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0}}
        wig = ["#comment",
               "track to skip",
               "chr1 0 5 1"
               ]
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")]}}
        all_gene_models = True
        self.assertEqual(coverageDict, readWiggle(wig, geneDict, all_gene_models))

    def testReadWiggleMulti(self):
        coverageDict = {"chr1": {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0, 15: 0, 16: 0, 17: 0, 18: 0,
                                 19: 0},
                        "chr2": {20: 0, 21: 7, 22: 7, 23: 7, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        wig = ["#comment",
               "track to skip",
               "chr1 0 5 1",
               "chr2\t21\t24 7",
               "chr3 0 10 2"
               ]
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}
                    }
        all_gene_models = True
        self.assertEqual(coverageDict, readWiggle(wig, geneDict, all_gene_models))

    def testReadAnnotationEmpty(self):
        geneDict = {}
        gtf = []
        bioType = None
        geneType = None
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType))
        gtf = ["# comment",
               "chr1\t\tnot_exon"
               ]
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType)),
        bioType = "good_bio"
        gtf = ["# comment",
               "chr1\t\tnot_exon",
               "chr1\tbad_bio\texon"
               ]
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType))

    def testReadAnnotationSingle(self):
        geneDict = {"gene1": {"tr1": [("chr2", 0, 10, "+")]}}
        gtf = ['# comment',
               'chr1\t\tnot_exon',
               'chr2\t\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr1";'
               ]
        bioType = None
        geneType = None
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType))

    def testReadAnnotationMulti(self):
        geneDict = {"gene1": {"tr1": [("chr2", 0, 10, "+"),
                                      ("chr2", 20, 50, "+")],
                              "tr2": [("chr2", 0, 10, "+")]},
                    "gene2": {"tr1": [("chr3", 0, 80, "-")]}
                    }
        gtf = ['chr2\t\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t20\t50\t\t+\t\tgene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr2";',
               'chr3\t\texon\t0\t80\t\t-\t\tgene_id "gene2";transcript_id "tr1";'
               ]
        bioType = None
        geneType = None
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType))

    def testReadAnnotationBioType(self):
        geneDict = {"gene1": {"tr1": [("chr2", 0, 10, "+"),
                                      ("chr2", 20, 50, "+")]},
                    "gene2": {"tr1": [("chr3", 0, 80, "-")]}
                    }
        gtf = ['chr2\tgoodType\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr1";',
               'chr2\tgoodType\texon\t20\t50\t\t+\t\tgene_id "gene1";transcript_id "tr1";',
               'chr2\tbad_type\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr2";',
               'chr3\tgoodType\texon\t0\t80\t\t-\t\tgene_id "gene2";transcript_id "tr1";'
               ]
        bioType = "goodType"
        geneType = None
        self.assertEqual(geneDict, readAnnotation(gtf, bioType, geneType))

    def testGenerateOutputArrayEmpty(self):
        outputArray = numpy.zeros(shape=100)
        geneDict = {}
        coverageDict = {}
        minGeneLength = 0
        outputfilename = None
        maxGeneLength = None
        doPrintList = False
        generatedArray  = createCoverageArray(
            geneDict, coverageDict,
            minGeneLength, maxGeneLength,
            outputfilename,
            doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)

    def testGenerateOutputArrayLengthFails(self):
        outputArray = numpy.zeros(shape=100)
        coverageDict = {"chr1": {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0, 15: 0, 16: 0, 17: 0, 18: 0,
                                 19: 0},
                        "chr2": {20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}}
        minGeneLength = 100
        verbose = False
        outputfilename = None
        maxGeneLength = None
        doPrintList = False
        generatedArray = createCoverageArray(geneDict, coverageDict,
                                             minGeneLength, maxGeneLength,
                                             outputfilename,
                                             doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)
        minGeneLength = 0
        maxGeneLength = 5
        generatedArray = createCoverageArray(geneDict, coverageDict,
                                             minGeneLength, maxGeneLength,
                                             outputfilename,
                                             doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)



if __name__ == '__main__':
    main()
