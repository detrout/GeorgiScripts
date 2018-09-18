import unittest
from SAM_reads_in_genes3_BAM import *


class SAMReadsInGenesBAMTest(unittest.TestCase):

    def testGetChromInfoEmpty(self):
        chrominfo = []
        chromInfoList = []
        self.assertEquals(chromInfoList, getChromInfo(chrominfo))

    def testGetChromInfoSingle(self):
        chrominfo = ["chr1\t30"]
        chromInfoList = [("chr1", 0, 30)]
        self.assertEquals(chromInfoList, getChromInfo(chrominfo))

    def testGetChromInfoMulti(self):
        chrominfo = ["chr1\t30", "chr2\t10"]
        chromInfoList = [("chr1", 0, 30), ("chr2", 0, 10)]
        self.assertEquals(chromInfoList, getChromInfo(chrominfo))

    def testGetReadMultiplicityEmpty(self):
        pass

    def testPileUpReadsEmpty(self):
        pass

    def testGetCountsEmpty(self):
        exonicReads = 0
        intronicReads = 0
        posCountsDict = {}
        geneDict = {}
        self.assertEqual((exonicReads, intronicReads),
                         getCounts(posCountsDict, geneDict, verbose=False))

    def testGetCountsSingle(self):
        geneDict = {"chr1": {"gene1": set([(1, 5), (15, 30)])}}
        exonicReads = 1
        intronicReads = 0
        posCountsDict = {"chr1": {20: 1}}
        self.assertEqual((exonicReads, intronicReads),
                         getCounts(posCountsDict, geneDict, verbose=False))
        exonicReads = 0
        intronicReads = 1
        posCountsDict = {"chr1": {10: 1}}
        self.assertEqual((exonicReads, intronicReads),
                         getCounts(posCountsDict, geneDict, verbose=False))
        exonicReads = 0
        intronicReads = 0
        posCountsDict = {"chr1": {40: 1}}
        self.assertEqual((exonicReads, intronicReads),
                         getCounts(posCountsDict, geneDict, verbose=False))

    def testGetCountsMulti(self):
        geneDict = {"chr1": {"gene1": [(1, 5), (15, 30)]}}
        exonicReads = 3
        intronicReads = 4
        posCountsDict = {"chr1": {3: 2, 10: 4, 20: 1, 40: 7}}
        self.assertEqual((exonicReads, intronicReads),
                         getCounts(posCountsDict, geneDict, verbose=False))

    def testReadGtfEmpty(self):
        gtf = []
        self.assertEquals({}, readGtf(gtf))

    def testReadGtfMulti(self):
        gtf = ['chr2\t\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t20\t50\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr1";',
               'chr2\t\texon\t0\t10\t\t+\t\t'
               'gene_id "gene1";transcript_id "tr2";',
               'chr3\t\texon\t0\t80\t\t-\t\t'
               'gene_id "gene2";transcript_id "tr1";'
               ]
        geneDict = {'chr2': {'gene1': set([(0, 10), (20, 50)])},
                    'chr3': {'gene2': set([(0, 80)])}
                    }
        self.assertEquals(geneDict, readGtf(gtf))

    def testWriteOutput(self):
        distribution = pseudoFile([])
        totalReads = 10.
        exonicReads = 5.
        intronicReads = 3.
        writeOutput(distribution, totalReads, exonicReads, intronicReads)
        self.assertEquals(distribution.read(), "Intergenic:0.2\t\n")
        self.assertEquals(distribution.read(), "Intronic:0.3\t\n")
        self.assertEquals(distribution.read(), "Exonic:0.5\t\n")
        self.assertEquals(distribution.read(), "#Class\tFraction\n")
