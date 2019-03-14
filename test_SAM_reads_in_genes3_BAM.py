import unittest
import os
from subprocess import check_call
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from SAM_reads_in_genes3_BAM import *


class SAMReadsInGenesBAMTest(unittest.TestCase):

    def testGetChromInfoEmpty(self):
        chrominfo = []
        chromInfoList = []
        self.assertEqual(chromInfoList, getChromInfo(chrominfo))

    def testGetChromInfoSingle(self):
        chrominfo = ["chr1\t30"]
        chromInfoList = [("chr1", 0, 30)]
        self.assertEqual(chromInfoList, getChromInfo(chrominfo))

    def testGetChromInfoMulti(self):
        chrominfo = ["chr1\t30", "chr2\t10"]
        chromInfoList = [("chr1", 0, 30), ("chr2", 0, 10)]
        self.assertEqual(chromInfoList, getChromInfo(chrominfo))

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
        self.assertEqual({}, readGtf(gtf))


    def testReadGtfComment(self):
        gtf = [
            '# test',
            '# more tests',
        ]
        self.assertEqual({}, readGtf(gtf))

    def testReadGtfNoExons(self):
        gtf = [
            'chr2\t\tCDS\t0\t10\t\t+\t\tgene_id "gene1";transcript_id "tr1";',
        ]
        self.assertEqual({}, readGtf(gtf))

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
        self.assertEqual(geneDict, readGtf(gtf))

    def testWriteOutput(self):
        distribution = pseudoFile([])
        totalReads = 10.
        exonicReads = 5.
        intronicReads = 3.
        writeOutput(distribution, totalReads, exonicReads, intronicReads)
        self.assertEqual(distribution.read(), "Intergenic:\t0.2\n")
        self.assertEqual(distribution.read(), "Intronic:\t0.3\n")
        self.assertEqual(distribution.read(), "Exonic:\t0.5\n")
        self.assertEqual(distribution.read(), "#Class\tFraction\n")

    def test_get_read_multiplicity_unique(self):
        chromInfoList = [("chr1", 0, 3000)]
        with prepare_bam('unique-read.sam') as unique:
            alignment = pysam.Samfile(unique, "rb")

            self.assertEqual(getReadMultiplicity(chromInfoList, alignment),
                             {'r0001': 1})

    def test_get_read_multiplicity_multi(self):
        chromInfoList = [("chr1", 0, 3000)]
        with prepare_bam('multireads.sam') as multireads:
            alignment = pysam.Samfile(multireads, "rb")

            self.assertEqual(getReadMultiplicity(chromInfoList, alignment),
                             {'r0001': 3})

    def test_pileup_reads_unique(self):
        chromInfoList = [("chr1", 0, 3000)]
        multiplicity = {'r0001': 1}
        with prepare_bam('unique-read-nh.sam') as unique:
            alignment = pysam.Samfile(unique, "rb")

            counts, reads = pileUpReads(
                chromInfoList, alignment, multiplicity)
            self.assertEqual(counts, {'chr1': {0: 1}})
            self.assertEqual(reads, 1)

    def test_pileup_reads_mutli(self):
        chromInfoList = [("chr1", 0, 3000)]
        multiplicity = {'r0001': 3}
        with prepare_bam('multireads-nh.sam') as multireads:
            alignment = pysam.Samfile(multireads, "rb")

            counts, reads = pileUpReads(
                chromInfoList, alignment, multiplicity)
            self.assertEqual(counts, {
                'chr1': {
                    0: 0.3333333333333333,
                    10: 0.3333333333333333,
                    20: 0.3333333333333333}})
            self.assertEqual(reads, 1)

    def test_pileup_reads_mutli(self):
        chromInfoList = [("chr1", 0, 3000)]
        multiplicity = {'r0001': 3}
        with prepare_bam('multireads.sam') as multireads:
            alignment = pysam.Samfile(multireads, "rb")

            self.assertRaises(SystemExit,
                              pileUpReads,
                              chromInfoList,
                              alignment,
                              multiplicity)

@contextmanager
def prepare_bam(samfile):
    with TemporaryDirectory('_sam_reads') as tempdir:
        try:
            path, name = os.path.split(samfile)
            bamfile = os.path.join(tempdir, name.replace('.sam', '.bam'))
            baifile = bamfile + '.bai'

            if not os.path.exists(bamfile):
                check_call(['samtools', 'view', '-t', 'bam', '-o', bamfile, samfile])
            if not os.path.exists(baifile):
                check_call(['samtools', 'index', bamfile])

            yield bamfile
        finally:
            if os.path.exists(bamfile):
                os.unlink(bamfile)
            if os.path.exists(baifile):
                os.unlink(baifile)

