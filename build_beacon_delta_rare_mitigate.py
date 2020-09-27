##Editor: Diyue Bu####
###last edition: 09/29/2015###

from copy import deepcopy
from math import log10
from sys import exit
from os import walk
from multiprocessing import Pool, Process
from random import seed, shuffle, random
# from numpy import random
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def personlist():
    ##remove family numbers
    frelate = '../../1000genomes/20140625_related_individuals.txt'
    f = open(frelate, 'r')
    lremove = []
    for l in f:
        if 'Sample' in l:
            continue
        lremove.append(l.strip().split()[0])
    f.close()

    fpinfo = '../../1000genomes/integrated_call_samples.20130502.ALL.ped'
    lpersonid = []
    lrelate = set([])
    f = open(fpinfo, 'r')
    for l in f:
        if 'Population' in l:
            continue
        r = l.strip().split()
        pid = r[1]
        if pid in lremove or pid in lrelate:
            continue
        relation = r[7]
        if relation == 'child':
            pass
        elif relation == 'mother' or relation == 'father':
            if r[8] == '0':
                lpersonid.append(pid)
            else:
                pass
        else:
            lpersonid.append(pid)
    del lrelate
    f.close()

    return lpersonid


def readcol(sfile):
    popSize = int(sfile.split('/')[-2].split('Size')[1])

    loutput = []

    with open(sfile, 'r') as f:
        SNP_count=0
        for l in f:
            l = l.strip()
            if l.startswith('##'):
                continue

            if l.startswith('#CHROM'):
                continue

            r = l.split()

            linfo = r[7].split(';')

            for i, elem in enumerate(linfo):
                if 'AC=' in elem:
                    index = i

            saf = linfo[index + 1].split('=')[1]
            if ',' in saf:
                laf = saf.split(',')
                faf = float(laf[0])
            else:
                faf = float(saf)
            # print faf

            # lsnpaf.append(faf)


            c = 0
            # check total # of snps in test set
            m = 0
            faf_release=faf
            for i in range(popSize):
                r_allele=r[9 + i].split('|')
                if '1' in r_allele[0]:
                    m += 1
                if '1' in r_allele[1]:
                    m += 1
            if m < t:
                x = 0
                faf_release = 0.0
            elif m == t:
                freq = random()
                if freq < delta:
                    x = 0
                    faf_release = 0.0
                else:
                    x = 1
            else:
                x = 1

            ######## Precomputing the Loglikelihood ratio ####################
            psuedo = delta #float(pow(10, -290)) #float(pow(10, -6))
            N = popSize #size of population
            D1 = dprob1(faf, t, N, delta)
            D0 = dprob0(faf, t, N)
            if D1 == 0.0:
                LH1 = x * log10(1 - D1 - psuedo) + (1 - x) * log10(D1 + psuedo)
            elif D1 == 1.0:
                LH1 = x * log10(1 - D1 + psuedo) + (1 - x) * log10(D1 - psuedo)
            else:
                LH1 = x * log10(1 - D1) + (1 - x) * log10(D1)
            if D0 == 0.0:
                LH0 = x * log10(1 - D0 - psuedo) + (1 - x) * log10(D0 + psuedo)
            elif D0 == 1.0:
                LH0 = x * log10(1 - D0 + psuedo) + (1 - x) * log10(D0 - psuedo)
            else:
                LH0 = x * log10(1 - D0) + (1 - x) * log10(D0)
            L = LH0 - LH1
            SNP_count += 1
            fpool = float(m)/2/float(N)
            re_diff = fpool-faf
            abs_diff = abs(re_diff)
            loutput.append([x, faf_release, faf, abs_diff, re_diff, L, (-1)*SNP_count])
    print '#SNP=' + str(SNP_count)
    return loutput


def main():
    ###walk file names and pair every two files to multi-process list####
    lfName = []
    for root, dirs, files in walk('../1000genomes/chr_vcf/'):
        for fname in files:
            lfName.append(fname)

    return lfName

##################################################
############# Added functions ####################
##################################################
def dprob1(f, t, N, delta):
    if t == 1:
        D = delta * pow((1 - f), (2 * (N - 1)))
        return D

    D = delta * pow((1 - f), (2 * (N - t)))
    for i in range(1, t):
        D += pow(f, (i - 1)) * pow((1 - f), (2 * (N - i)))

    return D

def dprob0(f, t, N):
    D = 0
    for i in range(t):
        D += pow(f, i) * pow((1 - f), (2 * (N - i)))

    return D


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-t', '--threshold')
    parser.add_argument('-d', '--delta')
    args = parser.parse_args()

    t = int(args.threshold)
    delta = float(args.delta)

    l10 = readcol('data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf')

    fname = 'data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf'
    schr = fname.split('/')[-1].split('chr')[1].split('_')[0]
    fout = open('chr' + schr + 'returnValue_500wholeR.txt', 'w')
    fout.write('returnValue\tminorAlleleFreq\n')
    ####################################
    ####### Mitigation Algorithm #######
    ####################################

    #### Configuration ######
    #### 1. Sorting on what ####
    Sort_base = 4 #1:allele freq. 2:The absolute difference between the pool and the underlying reference. 3: The difference between the pool and the underlying reference. 4:Loglikelihood ratio. 5:(-1) * default order
    #### 2. Configure Mitigation Strategy #####
    Mitigation_strategy = 1 #1:Only flip top K SNPs. 2. Only flip top SNPs smaller than a threshold Theta. 3: Only flip top k% SNPs.4:Automatically determined
    K = 100000
    k = 0.0001
    Theta = 0.05
    #### 3. Configure Noise Adding ########
    Noise_adding = 0 #1. Adding noise. 2. No noise.
    Noise_level = 0.001
    m_SNP = 0
    List_item = [];
    for lpos in l10:
        m_SNP += 1
        List_item.append(lpos[Sort_base+1])
    Sorted_List_item = sorted(List_item, reverse=True)
    if Mitigation_strategy == 1:
        Theta = Sorted_List_item[K]
    elif Mitigation_strategy == 3:
        Theta = Sorted_List_item[int(m_SNP*k)]
    elif Mitigation_strategy == 4:
        kappa=0.5
        fin = open('avg_powers.txt', 'r')
        for l in fin:
            Theta = Sorted_List_item[int(m_SNP*kappa*float(l.split('\n')[0]))]
    for lpos in l10:
        if lpos[0] == 0:
            print lpos
            if (Noise_adding == 1) & (random() < Noise_level):
                lpos[0] = 1
                lpos[1] = lpos[2]
        else:
            if lpos[Sort_base+1] >= Theta:
                if (Noise_adding == 1) & (random() < Noise_level):
                    pass
                else:
                    lpos[0] = 0
                    lpos[1] = 0
        fout.write(str(lpos[0]) + '\t' + str(lpos[1]) + '\n')
    fout.close()
