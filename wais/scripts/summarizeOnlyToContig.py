#!/usr/bin/python3

import argparse
import re
import sys

Col_qSeqId = 0
Col_sSeqId = 1
Col_pIdent = 2
Col_alignLen = 3
Col_sAlignStart = 8
Col_sAlignEnd = 9
Col_qAlignStart = 6
Col_qAlignEnd = 7
Col_qSeqLen = 12
Col_sSeqLen = 13

################################################## TOP_LVL
def summarizeTheRes(fn_only, fn_flanks, th_minPident, th_alignAndLenDiff):

    dict_only = loadTheOnly(fn_only, th_minPident, th_alignAndLenDiff)

    checkWithTheFlanks(fn_flanks, dict_only, th_minPident, th_alignAndLenDiff)

    #for readId in dict_only:
    #   print (readId + "\t" + str(dict_only[readId]))



################################################## AUX
def checkWithTheFlanks(fn_flanks, dict_only, th_minPident, th_alignAndLenDiff):
    print("ReadId\tContigInFlank(alignLen)\tisContigInPair\tContigInPair(alignLen)")
    with open(fn_flanks, 'r') as fh:
        for line in fh:
            if not re.match('^\#', line):
                line = line.strip()

                arr = line.split("\t")

                arr_readId = arr[Col_qSeqId].split("__")
                # print (arr_readId[0])
                if float(arr[Col_pIdent]) >= th_minPident:
                    if abs(int(arr[Col_alignLen]) - int(arr[Col_qSeqLen])) < th_alignAndLenDiff:

                        sys.stdout.write(arr[Col_qSeqId] + "\t" + arr[Col_sSeqId] + '(' + arr[Col_alignLen] + ')\t')

                        if (arr_readId[0] in dict_only):
                            # print (str(dict_only[arr_readId[0]]) + "\t" + line)
                            isContigPres = isContigPresent(dict_only[arr_readId[0]], arr[Col_sSeqId])


                            sys.stdout.write(str(isContigPres) + "\t")

                            #if isContigPres:

                            for (contigId_complRead, alignStartContig, alignEndContig) in dict_only[arr_readId[0]]:
                                alignLen = max([int(alignEndContig),int(alignStartContig)]) - min([int(alignEndContig), int(alignStartContig)]) + 1

                                sys.stdout.write(contigId_complRead + "(" + str(alignLen) + ");")



                        sys.stdout.write('\n')

                            # print("\tFound")


def isContigPresent(arr_contigAligns, contigName):

    for (contigId_complRead, alignStartContig, alignEndContig) in arr_contigAligns:
        if contigName == contigId_complRead:
            return True

    return False


def loadTheOnly(fn_only, th_minPident, th_alignAndLenDiff):

    dict_only = dict() # dict_{readId} => [(contigName, alignStart, alignEnd), (contigName, alignStart, alignEnd), ...]

    with open(fn_only, 'r') as fh:

        for line in fh:
            if not re.match('^\#', line):
                arr = line.split("\t")

                if float(arr[Col_pIdent]) >= th_minPident:
                    if abs(int(arr[Col_alignLen]) - int(arr[Col_qSeqLen])) < th_alignAndLenDiff:
                        # print (arr[Col_qSeqId] + "\t" + arr[Col_sSeqId] + "\t" + arr[Col_alignLen] + "\t" + arr[Col_qSeqLen])

                        # Add to dict
                        if arr[Col_qSeqId] not in dict_only:
                            dict_only[arr[Col_qSeqId]] = []

                        dict_only[arr[Col_qSeqId]].append((arr[Col_sSeqId], arr[Col_sAlignStart], arr[Col_sAlignEnd]))

    return dict_only


################################################## MAIN
def main():
    parser = argparse.ArgumentParser(description='Check the alignments of only reads to contigs, that its between 300bps of its pairs')

    parser.add_argument('--only', nargs=1, required=True, help="Blast results of only1|2.fasta to contigs")
    parser.add_argument('--flanks', nargs=1, required=True, help="Blast results of partial aligned reads (aka. flanks, aka pairs) to contigs")

    parser.add_argument('--th_minPident', nargs=1, type=float, default=[0], help="Minimum percent identity of alignment to keep (flanks to contig)")
    parser.add_argument('--th_alignAndLenDiff', nargs=1, type=int, default=[10000], help="Alignment and length difference to keep (flanks to contig)")


    args = parser.parse_args()

    summarizeTheRes(args.only[0], args.flanks[0], args.th_minPident[0], args.th_alignAndLenDiff[0])

if __name__ == '__main__':
    main()
