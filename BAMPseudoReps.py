##################################
#                                #
# Last modified 03/09/2015       #
#                                #
# Georgi Marinov                 #
#                                #
##################################
from __future__ import print_function

import sys
import os
import random
import pysam
import string


def main():
    if len(sys.argv) < 2:
        print('usage: python %s BAM' % sys.argv[0])
        print('Only run this script for files containing uniquely aligned'
              'reads; the script will not check for alignment multiplicity!')
        sys.exit(1)

    BAM = sys.argv[1]
    BAM_basename = os.path.splitext(BAM)[0]

    samfile = pysam.Samfile(BAM, "rb")

    outfile1 = pysam.Samfile(BAM_basename + '.pseudoRep1.bam', "wb",
                             template=samfile)
    outfile2 = pysam.Samfile(BAM_basename + '.pseudoRep2.bam', "wb",
                             template=samfile)

    for alignedread in samfile.fetch():
        chr = samfile.getrname(alignedread.tid)
        if chr == '*':
            continue
        if random.random() >= 0.5:
            outfile1.write(alignedread)
        else:
            outfile2.write(alignedread)

    outfile1.close()
    outfile2.close()

if __name__ == '__main__':
    main()
