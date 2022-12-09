#! /usr/bin/python3

import sys
import math
import argparse
from collections import Counter
from scipy.spatial import distance
from parse_fasta import parse_fasta


def get_columns(seqDict):
    colDict = {}
    pos = 0
    for k, v in seqDict.items():
        for i in v:
            try:
                colDict[pos].append(i)
            except KeyError:
                colDict[pos] = []
                colDict[pos].append(i)
            pos += 1
        pos = 0
    return colDict


def calc_col_prop(colDict):
    colPropDict = {}
    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V"]
    for k, v in colDict.items():
        tot = sum([x[1] for x in Counter(v).items() if x[0] != "-"])
        props = []
        for i in aa:
            try:
                c = Counter(v)[i]
            except KeyError:
                c = 0
            props.append(c / tot)
            colPropDict[k] = props
    return colPropDict


def calc_euclidean(colPropDict1, colPropDict2):
    distDict = {}
    for k, v in colPropDict1.items():
        distDict[k] = sum([(a - b)**2 for a, b in zip(v,
                                                      colPropDict2[k])])**0.5
    return distDict


def calc_kl(colProp1, colProp2):
    div = [x * math.log2((x / y)) for x, y in zip(colProp1, colProp2)
           if y != 0 and x != 0]
    kl = sum(div)
    return kl


def calc_jsd(colPropDict1, colPropDict2):
    distDict = {}
    for k, v in colPropDict1.items():
        m = [(x + y) / 2 for x, y in zip(v, colPropDict2[k])]
        jsd = (calc_kl(v, m) / 2) + (calc_kl(colPropDict2[k], m) / 2)
        distDict[k] = jsd
    return distDict


def calc_jsd_scipy(colPropDict1, colPropDict2, base=2):
    distDict = {}
    for k, v in colPropDict1.items():
        distDict[k] = distance.jensenshannon(v, colPropDict2[k],
                                             base=base) ** 2
    return distDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", help="Alignment in FASTA format")
    parser.add_argument("flag1", help="Flag to recognise sequences \
                        in group 1. Must be first substring of name.")
    parser.add_argument("flag2", help="Flag to recognise sequences \
                        in group 2. Must be first substring of name.")
    parser.add_argument("-d", "--distance", help="distance \
                        metric to use [euc (default)/jsd]",
                        default="euc")
    args = parser.parse_args()

    seqs1 = dict([x for x in parse_fasta(args.sequence)
                 if x[0].startswith(args.flag1)])
    # for k in seqs1.keys():
    #     print(k)
    seqs2 = dict([x for x in parse_fasta(args.sequence)
                 if x[0].startswith(args.flag2)])
    # for k in seqs2.keys():
    #     print(k)
    cols1 = get_columns(seqs1)
    cols2 = get_columns(seqs2)
    colProps1 = calc_col_prop(cols1)
    # print(colProp1)
    colProps2 = calc_col_prop(cols2)
    # print(colProp2)
    if args.distance == "euc":
        dist = calc_euclidean(colProps1, colProps2)
    elif args.distance == "jsd":
        dist = calc_jsd_scipy(colProps1, colProps2)
    else:
        print("distance argument not recognised")
        sys.exit()
    print("pos,dist")
    tmp = [(a+1, b) for a, b in zip(sorted(dist, key=dist.get, reverse=True),
           sorted(dist.values(), reverse=True))]
    for i in tmp:
        print(str(i[0])+","+str(i[1]))
