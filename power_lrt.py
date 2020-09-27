##########################USAGE##################################################
# this script is to calculate attack power based on lrt scores			#
# python figure_beacondb_retry.py						#
# input: lrt scores from previous step (need to specify input files in script)	#
# output: figure w/ specified inputs (e.g. power of different error rates)	#
#										#
#################################################################################

import os
import matplotlib.pyplot as plt
from random import sample
import sys
from copy import deepcopy
from multiprocessing import Pool, Process
import numpy


def substract(l_co, l_test):
    # print len(l_co), len(l_test)
    dPowerThr = 0.95
    # nmaxSize = len(l_co)
    # L = []
    nthr = sorted(l_co)[int(len(l_co) * (1 - dPowerThr))]
    # print 'nthr', nthr
    ncase = 0
    for j in sorted(l_test):
        #	print j
        if j < nthr:
            ncase += 1

    power = float(ncase) / float(len(l_test))

    # print lper
    return power


def median(ltarget):
    # print type(ltarget)
    lsort = sorted(ltarget)
    nl = len(lsort)
    if nl % 2 == 0:
        nmedian = lsort[nl / 2]
    else:
        ni = int(float(nl / 2) + 1)
        nmedian = lsort[ni]
    return nmedian


def main1(delta):
    Lcase = []
    Ltest = []
    b = 0
    fcase = 0
    ftest = 0
    # for root, dirs, files in os.walk('results_beacondb/t' + str(t) + '/'):
    #        for fname in files:
    fname = 'chr10_lrtScoreByStep.txt'
    cmin = 100000
    c = 0
    with open('results_beacondb/Raredelta_largef' + str(delta) + '/' + fname, 'r') as f:

        #                f = open('results_beacondb/t' + str(t) + '/' + fname, 'r')
        for l in f:
            if 'query' in l or 'Not' in l:
                continue

            if 'individuals' in l:
                if c != 0 and c < cmin:
                    cmin = c
                c = 0
                continue
            c += 1
    if c < cmin:
        cmin = c
    f.close()

    Lcase = [[] for i in range(cmin)]
    Ltest = [[] for i in range(cmin)]

    with open('results_beacondb/Raredelta_largef' + str(delta) + '/' + fname, 'r') as f:
        for l in f:
            if 'query' in l or 'Not' in l:
                continue

            if 'individuals' in l:
                fcase = 0
                ftest = 0
                c = 0
                continue
            r = l.strip().split()
            fcase += float(r[0])
            ftest += float(r[1])
            c += 1
            try:
                Lcase[int(c - 1)].append(fcase)
                Ltest[int(c - 1)].append(ftest)
            except:
                pass

    f.close()
    lpower = []
    print cmin
    for i in range(cmin):
        npower = substract(Lcase[i], Ltest[i])
        lpower.append(npower)
    return lpower


def main(delta):
    Lcase = []
    Ltest = []
    b = 0
    fcase = 0
    ftest = 0
    #############################################
    #  lrt score output file, specified by user #
    #############################################
    fname = 'chr10_lrtScoreByStep.txt'
    cmin = 100000
    c = 0
    with open(fname, 'r') as f:
        for l in f:
            if 'query' in l or 'Not' in l:
                continue

            if 'individuals' in l:
                if c != 0 and c < cmin:
                    cmin = c
                c = 0
                continue
            c += 1
    if c < cmin:
        cmin = c
    f.close()

    Lcase = [[] for i in range(cmin)]
    Ltest = [[] for i in range(cmin)]

    with open(fname, 'r') as f:
        for l in f:
            if 'query' in l or 'Not' in l:
                continue

            if 'individuals' in l:
                fcase = 0
                ftest = 0
                c = 0
                continue
            r = l.strip().split()
            fcase += float(r[0])
            ftest += float(r[1])
            c += 1
            try:
                Lcase[int(c - 1)].append(fcase)
                Ltest[int(c - 1)].append(ftest)
            except:
                pass

    f.close()

    lpower = []
    b = False
    print cmin
    for i in range(cmin):
        npower = substract(Lcase[i], Ltest[i])
        # print npower
        if npower == 1.0 and b == False:
            # print i
            b = True
        else:
            pass
        lpower.append(npower)
    return lpower


def draw_pernum(Lpower, lper):
    fout = open('avg_powers.txt', 'w') #log result
    # print lper
    lcolor = ['xb-', 'xg-', 'xr-', 'xc-', 'xm-', 'xy-', 'xk-', 'x#b3de69-', 'x#A60628-', 'x#988ED5-', '#7A68A6--']
    plt.figure()
    plt.hold(True)
    ncounter = 0
    for i in range(len(Lpower)):
        powers = Lpower[i]
        print powers
        fout.write(str(numpy.mean(powers)) + '\n') #log result
        delta = ldelta[i]
        ncounter += 1
        plt.plot(lper, powers, color=lcolor[ncounter - 1][1:-1], ls='-', label='delta=' + str(delta))
    fout.close()
    l5 = [0.05 for i in range(len(lper))]
    plt.plot(lper, l5, color='r', ls='-.', label='5% false positive rate line')

    plt.ylim([0, 1.1])
    plt.xlabel('# of snps queried by user', fontsize=12)
    plt.ylabel('True positive of Likelihood Ratio Test')
    plt.legend(loc='center right', prop={'size': 10})
    plt.title('query # detect difference with different error rate on rare snps(sorted query)')

    ###################################
    #   figure name, specified by user#
    ###################################
    fname = 'power_delta.pdf'
    plt.savefig(fname, format='pdf')


if __name__ == '__main__':

    #############################################
    #   different error rate specified by user  #
    #############################################
    # ldelta = [1e-06, 0.0001, 0.001, 0.01, 0.05, 0.15, 0.5]
    ldelta = [1e-06]
    #####################################
    #  beacondb size, specified by user #
    #####################################
    N = 500
    # p = Pool(4)
    # Lpower = p.map(main,ldelta)
    Lpower = [main(ldelta[0])]

    length = 10000000
    for lpower in Lpower:
        if length > len(lpower):
            length = len(lpower)
    lper = [(i) * 100 + 1 for i in range(length)]
    for i in range(len(Lpower)):
        Lpower[i] = Lpower[i][:length]

    draw_pernum(Lpower, lper)
