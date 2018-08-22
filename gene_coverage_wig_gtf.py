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
import os
import unittest

logger = logging.getLogger('gene_coverage_wig_gtf')

try:
    import pyBigWig
except ImportError:
    logger.warn('pyBigWig is not available, no direct bigwig support')
    pyBigWig = None

version = '2.0'

NORMALIZATIONS = {
    'sum': numpy.mean,
    'max': numpy.max,
    'mean': numpy.mean,
    'median': numpy.median,
}

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.quiet:
        logging.basicConfig(level=logging.WARN)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.version:
        parser.exit(message="version {}\n".format(version))

    if args.source_type is not None:
        logger.info('will only consider genes from source %s', args.source_type)

    if args.max_gene_length is not None:
        logger.info('will only consider genes longer than %s and shorter than %s',
                    args.min_gene_length,
                    args.max_gene_length)
    else:
        logger.info('will only consider genes longer than %s', args.min_gene_length)

    if args.normalize:
        logger.info('will normalize scores')

    if args.gene_normalization in NORMALIZATIONS:
        logger.info('Will normalize gene bins by %s', args.gene_normalization)

    if args.all_gene_models:
        logger.info('will only all genes')
    else:
        logger.info('will only use genes with one isoform')

    if args.gene_type:
        logger.info('will only consider genes of type %s', args.gene_type)

    geneDict = loadAnnotation(args.gtf, args.source_type, args.gene_type)
    logger.info('genes passed type filters %s', len(geneDict))

    gene_coverage = loadGeneCoverage(args.filename, geneDict, args.all_gene_models)

    if args.print_list:
        geneListFilename = args.output + '.geneList'
    else:
        geneListFilename = None

    outputArray = createCoveragePercentiles(
        geneDict, gene_coverage,
        args.min_gene_length, args.max_gene_length,
        geneListFilename,
        args.normalize,
        args.gene_normalization,
    )

    with open(args.output, 'wt') as outfile:
        for i in range(100):
            outfile.write(str(i) + '\t' + str(outputArray[i])+'\n')

def make_parser():
    parser = ArgumentParser()
    parser.add_argument('--gtf', required=True, help='GTF file name')
    parser.add_argument('filename', help='file of type gedgraph or bigwig to compute coverage of')
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
                        help='normalize scores by number of genes')
    parser.add_argument('--gene-normalization',
                        choices=['none'] + sorted(NORMALIZATIONS.keys()),
                        default='none',
                        help='Per gene model normalization. sum is total number of reads assigned to gene. '\
                        'max is the maximum gene bin size.')

    parser.add_argument('-q', '--quiet', action='store_true',
                        help='only report errors')
    parser.add_argument('--version', action='store_true', help='report version number')
    return parser

def loadAnnotation(filename, source_type, gene_type_filter):
    with open(filename, 'rt') as gtfStream:
        return parseAnnotation(
            gtfStream,
            source_type,
            gene_type_filter)

def parseAnnotation(stream, source_type, gene_type_filter):
    geneDict={}
    entriesToDelete=set()
    for line in stream:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if source_type is not None and fields[1] != source_type:
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
    start = field.find(name)
    if start == -1:
        return None
    start += len(name)
    start = field.find('"', start)
    if start == -1:
        return None
    start += 1
    end = field.index('"', start)
    return field[start:end]


def loadGeneCoverage(filename, geneDict, all_gene_models):
    """Read coverage data for genes of interest

    :Parameters:
      - filename: source file to read
      - geneDict: GTF annotations
      - all_gene_models: flag if we should use all gene models instead
        of just single gene models
    """
    instream = guessFileOpen(filename)
    if instream is None:
        logger.error('Unable to open %s. Not a supported file-type', os.path.abspath(filename))
        raise RuntimeError('Unsupported file type')

    if pyBigWig and isinstance(instream, pyBigWig.pyBigWig) and instream.isBigWig():
        gene_coverage = readBigwig(instream, geneDict, all_gene_models)
    else:
        gene_coverage = readWiggle(instream, geneDict, all_gene_models)

    return gene_coverage

def guessFileOpen(filename):
    """Try opening the file with supported file readers

    returns the opened object
    """
    try:
        stream = open(filename, 'rt')
        header = stream.read(100)
        for c in header:
            if not (c.isprintable() or c.isspace()):
                logger.error('Not a text file %s %s', c, ord(c))
                raise ValueError('Not a text file %s %s', c, ord(c))
        stream.seek(0)
        return stream
    except ValueError as e:
        logger.warning('%s is not a text file. %s', os.path.abspath(filename), str(e))

    if pyBigWig:
        try:
            stream = pyBigWig.open(filename)
            return stream
        except RuntimeError as e:
            logger.debug('%s is not a bam file. %s', filename, str(e))

    if not os.path.exists(filename):
        logger.error('%s is not a local file', filename)

    return None


def readBigwig(bigwig, geneDict, all_gene_models):
    coverageDict = initializeCoverageDict(geneDict, all_gene_models)
    chromosomes = set(coverageDict).intersection(set(bigwig.chroms()))
    for chromosome in chromosomes:
        for left, right, score in bigwig.intervals(chromosome):
            for j in range(left, right):
                if j in coverageDict[chromosome]:
                    coverageDict[chromosome][j]=score

    logger.info('finished inputting bigwig')
    logger.info('genes passed type filters %s', len(geneDict))
    return coverageDict


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
    return coverageDict


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


def createCoveragePercentiles(GeneDict, coverageDict,
                        minGeneLength, maxGeneLength=None,
                        geneListFilename=None,
                        doNormalize=False,
                        gene_normalization=None):
    outputArray = numpy.zeros(shape=100)

    geneListStream = open(geneListFilename, 'wt') if geneListFilename else None

    geneNumber=0
    for geneID in GeneDict:
        NucleotideList, chromosome = buildNucleotideList(GeneDict[geneID])
        geneLength=len(NucleotideList)
        if geneLength < minGeneLength:
            continue
        if maxGeneLength is not None and geneLength > maxGeneLength:
            continue
        final_vector = numpy.zeros(shape=100)
        bins = numpy.linspace(0, geneLength, num = 101, dtype=int)
        start = bins[0]
        for i, end in enumerate(bins[1:]):
            counts = [coverageDict[chromosome][pos] for pos in NucleotideList[start:end]]
            final_vector[i] = numpy.mean(counts)
            start = end
        assert len(final_vector) == 100
        final_vector_sum = numpy.sum(final_vector)
        if final_vector_sum > 0:
            if gene_normalization in NORMALIZATIONS:
                final_vector /= NORMALIZATIONS[gene_normalization](final_vector)

            outputArray += final_vector
            geneNumber+=1

            if geneListStream:
                geneListStream.write(geneID)
                geneListStream.write('\t')
                geneListStream.write('\t'.join([str(x) for x in final_vector]))
                geneListStream.write('\n')

    logger.info('%s genes considered', geneNumber)

    if doNormalize:
        outputArray /= float(geneNumber)

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

    def testInitializeCoverageDictMultiAllModels(self):
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

    def testInitializeCoverageDictMultiSingleModel(self):
        coverageDict = {"chr2": {20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}}
        all_gene_models = False
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

    def testReadWiggleMultiAllModels(self):
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

    def testReadWiggleMultiSingleModels(self):
        coverageDict = {"chr2": {20: 0, 21: 7, 22: 7, 23: 7, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        wig = ["#comment",
               "track to skip",
               "chr1 0 5 1",
               "chr2\t21\t24 7",
               "chr3 0 10 2"
               ]
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}}
        all_gene_models = False
        self.assertEqual(coverageDict, readWiggle(wig, geneDict, all_gene_models))

    def testReadAnnotationEmpty(self):
        geneDict = {}
        gtf = []
        bioType = None
        geneType = None
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))
        gtf = ["# comment",
               "chr1\t\tnot_exon"
               ]
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType)),
        bioType = "good_bio"
        gtf = ["# comment",
               "chr1\t\tnot_exon",
               "chr1\tbad_bio\texon"
               ]
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))

    def testReadAnnotationSingle(self):
        geneDict = {"gene1": {"tr1": [("chr2", 0, 10, "+")]}}
        gtf = ['# comment',
               'chr1\t\tnot_exon',
               'chr2\t\texon\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr1";'
               ]
        bioType = None
        geneType = None
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))

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
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))

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
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))

    def testGenerateOutputArrayEmpty(self):
        outputArray = numpy.zeros(shape=100)
        geneDict = {}
        coverageDict = {}
        minGeneLength = 0
        outputfilename = None
        maxGeneLength = None
        doPrintList = False
        generatedArray  = createCoveragePercentiles(
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
        generatedArray = createCoveragePercentiles(geneDict, coverageDict,
                                             minGeneLength, maxGeneLength,
                                             outputfilename,
                                             doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)
        minGeneLength = 0
        maxGeneLength = 5
        generatedArray = createCoveragePercentiles(geneDict, coverageDict,
                                             minGeneLength, maxGeneLength,
                                             outputfilename,
                                             doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)


    def testComputeCoverageWithIntrons(self):
        outputArray = numpy.asarray([8 for x in range(100)])
        outputArray[50:60] = 4
        coverageDict = {"chr1": {}}
        # add exon coverage for two equally expressed transcripts, one
        # of which has an extra exon.
        coverageDict["chr1"].update({x: 8 for x in range(0, 500)})
        coverageDict["chr1"].update({x: 4 for x in range(600, 700)})
        coverageDict["chr1"].update({x: 8 for x in range(1600, 2000)})

        # add an express intron
        coverageDict["chr1"].update({x: 10 for x in range(900, 1100)})

        geneDict = {"gene1":
                    {"transcript1": [
                        ("chr1", 0, 500, "+"),
                        ("chr1", 1600, 2000, "+"),
                    ],
                     "transcript2": [
                         ("chr1", 0, 500, "+"),
                         ("chr1", 600, 700, "+"),
                         ("chr1", 1600, 2000, "+"),
                    ]}}
        minGeneLength = 100
        verbose = False
        outputfilename = None
        maxGeneLength = None
        doPrintList = False
        generatedArray = createCoveragePercentiles(
            geneDict, coverageDict,
            minGeneLength, maxGeneLength,
            outputfilename,
            doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)

if __name__ == '__main__':
    main()
