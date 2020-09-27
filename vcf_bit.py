#!/usr/bin/env python


import bitarray as ba
import pickle as pl
from multiprocessing import Pool, Process
from argparse import ArgumentParser


def main(sIn):
    lid = []
    ######################################################
    # change the number of individuals in input vcf file #
    ######################################################
    lRes = [ba.bitarray() for i in range(N)]
    count = 0
    with open(sIn, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue

            count += 1
            r = l.split()
            index = 0
            for gt in r[9:]:
                if '1' in gt:
                    lRes[index].append(True)
                else:
                    lRes[index].append(False)
                index += 1

    fname = sIn.split('/')[-2]
    output = open(fname + 'wholeRecords.pkl', 'wb')
    pl.dump(lRes, output, -1)
    output.close()


if __name__ == '__main__':
    ##########################
    #  change input vcf file #
    ##########################
    # fname = '/data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf'
    parser = ArgumentParser()
    parser.add_argument('-f', '--filename')
    parser.add_argument('-n', '--N')
    args = parser.parse_args()

    fname = args.filename
    N = int(args.N)

    main(fname)
