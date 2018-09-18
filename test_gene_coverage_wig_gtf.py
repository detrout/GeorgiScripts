import unittest
from gene_coverage_wig_gtf import *


class GeneCoverageWigGtfTest(unittest.TestCase):
    def testInitializeCoverageDictEmpty(self):
        coverageDict = {}
        geneDict = {}
        all_gene_models = True
        self.assertEqual(coverageDict,
                         initializeCoverageDict(geneDict, all_gene_models))

    def testInitializeCoverageDictSimple(self):
        coverageDict = {"chr1": {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")]}}
        all_gene_models = True
        self.assertEqual(coverageDict,
                         initializeCoverageDict(geneDict, all_gene_models))

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
        self.assertEqual(coverageDict,
                         initializeCoverageDict(geneDict, all_gene_models))

    def testInitializeCoverageDictMultiSingleModel(self):
        coverageDict = {"chr2": {20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0,
                                 26: 0, 27: 0, 28: 0, 29: 0}}
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")],
                              "transcript2": [("chr1", 15, 20, "+")]},
                    "gene2": {"transcript1": [("chr2", 20, 30, "+")]}}
        all_gene_models = False
        self.assertEqual(coverageDict,
                         initializeCoverageDict(geneDict, all_gene_models))

    def testReadWiggleEmpty(self):
        coverageDict = {}
        wig = []
        geneDict = {}
        all_gene_models = True
        self.assertEqual(coverageDict,
                         readWiggle(wig, geneDict, all_gene_models))

    def testReadWiggleSingle(self):
        coverageDict = {"chr1": {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 0, 6: 0,
                                 7: 0, 8: 0, 9: 0}}
        wig = ["#comment",
               "track to skip",
               "chr1 0 5 1"
               ]
        geneDict = {"gene1": {"transcript1": [("chr1", 0, 10, "+")]}}
        all_gene_models = True
        self.assertEqual(coverageDict,
                         readWiggle(wig, geneDict, all_gene_models))

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
        self.assertEqual(coverageDict,
                         readWiggle(wig, geneDict, all_gene_models))

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
        self.assertEqual(coverageDict,
                         readWiggle(wig, geneDict, all_gene_models))

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
               'chr2\t\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";'
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
        gtf = ['chr2\t\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t20\t50\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr2";',
               'chr3\t\texon\t0\t80\t\t-\t\t'
               'gene_id "gene2";transcript_id "tr1";'
               ]
        bioType = None
        geneType = None
        self.assertEqual(geneDict, parseAnnotation(gtf, bioType, geneType))

    def testReadAnnotationBioType(self):
        geneDict = {"gene1": {"tr1": [("chr2", 0, 10, "+"),
                                      ("chr2", 20, 50, "+")]},
                    "gene2": {"tr1": [("chr3", 0, 80, "-")]}
                    }
        gtf = ['chr2\tgoodType\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\tgoodType\texon\t20\t50\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\tbad_type\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr2";',
               'chr3\tgoodType\texon\t0\t80\t\t-\t\t'
               'gene_id "gene2";transcript_id "tr1";'
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
        generatedArray = createCoveragePercentiles(
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
                                                   minGeneLength,
                                                   maxGeneLength,
                                                   outputfilename,
                                                   doPrintList)
        numpy.testing.assert_array_equal(outputArray, generatedArray)
        minGeneLength = 0
        maxGeneLength = 5
        generatedArray = createCoveragePercentiles(geneDict, coverageDict,
                                                   minGeneLength,
                                                   maxGeneLength,
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
