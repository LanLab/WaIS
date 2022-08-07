#!/usr/bin/python3

import argparse
import re
import sys

from Bio import SeqIO

import getFlankingSeqs as fns_flankingSeqs

################################# TOP_LVL
def rmTheFlanks(blastRes_1, blastRes_2, minChoppedLen, fn_flanks_1, fn_flanks_2, fn_out_flanks1, fn_out_flanks2, fn_flanks1AndOnly2, fn_flanks2AndOnly1): #, fn_flanksOut_1, fn_flanksOut_2): # flanks_1, flanks_2, fn_out_1, fn_out_2):

    #dict_reads1 = loadTheReadNamesAndLen(fn_reads1)
    #dict_reads2 = loadTheReadNamesAndLen(fn_reads1)

    #for key in dict_reads1:
    #    print (key + "\t" + str(dict_reads1[key]))

    # list_contigsComplAligned = list()


    dict_readsToIS_1 = fns_flankingSeqs.loadBlastRes(blastRes_1)
    dict_readsToIS_2 = fns_flankingSeqs.loadBlastRes(blastRes_2)

    # doTheOutput(dict_readsToIS_1, dict_readsToIS_2, minChoppedLen, 'both.txt', 'only1.txt', 'only2.txt')
    (list_readsToKeep_1, list_readsToKeep_2) = doTheOutput(dict_readsToIS_1, dict_readsToIS_2, minChoppedLen, 'both.txt', 'only1.txt', 'only2.txt')
    
    printTheFlanks_(list_readsToKeep_1, fn_flanks_1, fn_out_flanks1)
    printTheFlanks_(list_readsToKeep_2, fn_flanks_2, fn_out_flanks2)
    
    # rewriteTheFlanks();

    """
    loadBlastRes(blastRes_1, list_contigsComplAligned)
    loadBlastRes(blastRes_2, list_contigsComplAligned)

    loadTheReadNamesAndLen(flanks_1, fn_out_1, list_contigsComplAligned)
    loadTheReadNamesAndLen(flanks_2, fn_out_2, list_contigsComplAligned)
    """

################################# AUX
def printTheFlanks_(list_readsToKeep, fn_flanks, fn_out_flanks):


    fh_out = open(fn_out_flanks, 'w+')


    # dict_readLens = dict() # dict_{readName} => len


    with open(fn_flanks, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            arr = record.id.split("__")

            if (arr[0] in list_readsToKeep):
                fh_out.write(">" + record.id + "\n")
                fh_out.write(str(record.seq) + "\n")
            # else:
            #    print (arr[0])

    fh_out.close() 





def getTheISAndInfo(tup):
    # print(tup)
    (alignStart, alignEnd, readLen, pre_, post_, alignDir) = tup
    alignLen = max([alignStart, alignEnd]) - min([alignStart, alignEnd])

    choppedLen = readLen - alignLen;
    return (alignLen, readLen, choppedLen)


def doTheOutput(dict_readsToIS_1, dict_readsToIS_2, minChoppedLen, fn_both, fn_only1, fn_only2):

    list_readsToKeep_1 = list() # [readName, readName, ...]
    list_readsToKeep_2 = list()

    with open(fn_both, 'w+') as fh_:
        fh_.write("ContigName\tISname1:alignLen1:readLen1:choppedLen1:isKeep1\tSname2:alignLen2:readLen2:choppedLen2:isKeep2\tSingle:Multiple\tAllSame\tToKeepRow|AddedTo1|AddedTo2\n")
        for contigName in dict_readsToIS_1:
            if contigName in dict_readsToIS_2:

                arr_ISin1 = []
                arr_ISin2 = []

                isAnyRm_1 = False
                isAnyRm_2 = False


                fh_.write(contigName + "\t")

                for ISname in dict_readsToIS_1[contigName]:
                    

                    for anInfo in dict_readsToIS_1[contigName][ISname]: 
                        # print (anInfo)
                        (alignLen, readLen, choppedLen) = getTheISAndInfo(anInfo)
                        # print (str(alignLen) + "\t" + str(readLen) + "\t" + str(choppedLen))
                    
                        isKeep="Rm"

                        if choppedLen >= minChoppedLen:
                            isKeep = "Keep"
                        else:
                            isAnyRm_1 = True


                        fh_.write(ISname + ":" + str(alignLen) + ":" + str(readLen) +  ":" + str(choppedLen) + ":" + isKeep +   ';')

                        arr_ISin1.append(ISname)

                fh_.write('\t')


                for ISname in dict_readsToIS_2[contigName]:
                    
                    for anInfo in dict_readsToIS_2[contigName][ISname]:
                        (alignLen, readLen, choppedLen) = getTheISAndInfo(anInfo)

                        isKeep="Rm"

                        if choppedLen >= minChoppedLen:
                            isKeep = "Keep"
                        else:
                            isAnyRm_2 = True


                        fh_.write(ISname + ":" + str(alignLen) + ":" + str(readLen) +  ":" + str(choppedLen) + ":" + isKeep +  ';')

                        arr_ISin2.append(ISname)

                fh_.write('\t')


                isAnyNotSingle = False
                if len(dict_readsToIS_1[contigName]) > 1:
                    fh_.write('Multiple')
                    isAnyNotSingle = True
                else:
                    fh_.write("Single")

                if len(dict_readsToIS_2[contigName]) > 1:
                    fh_.write(':Multiple')
                    isAnyNotSingle = True
                else:
                    fh_.write(":Single")

                fh_.write('\t')


                areAllSame = False
                if len(arr_ISin1) == len(arr_ISin2):
                    intersecting = set(arr_ISin1).intersection(set(arr_ISin2))


                    if len(intersecting) == len(arr_ISin1):
                        fh_.write("AllSame")
                        areAllSame = True
                fh_.write('\t')

                if isAnyNotSingle == False and areAllSame == True:
                    fh_.write('ToKeepRow')

                    if isAnyRm_1 == False:
                        list_readsToKeep_1.append(contigName)
                        fh_.write('|1')
                    if isAnyRm_2 == False:
                        list_readsToKeep_2.append(contigName)
                        fh_.write('|2')



                fh_.write('\n')



    print("######### FOUND IN /1 only")
    with open(fn_only1, 'w+') as fh_:
        for contigName in dict_readsToIS_1:
            if contigName not in dict_readsToIS_2:
                fh_.write(contigName + '\t')
                for ISname in dict_readsToIS_1[contigName]:
                    
                    for anInfo in dict_readsToIS_1[contigName][ISname]:
                        (alignLen, readLen, choppedLen) = getTheISAndInfo(anInfo)

                        isKeep="Rm"

                        if choppedLen >= minChoppedLen:
                            isKeep = "Keep"
                        else:
                            isAnyRm_2 = True


                        fh_.write(ISname + ":" + str(alignLen) + ":" + str(readLen) +  ":" + str(choppedLen) + ":" + isKeep +  ';')

                        arr_ISin2.append(ISname)

                fh_.write('\t')


                isAnyNotSingle = False
                if len(dict_readsToIS_1[contigName]) > 1:
                    fh_.write('Multiple')
                    isAnyNotSingle = True
                else:
                    fh_.write("Single")

                fh_.write('\t')
                if isKeep == "Keep" and isAnyNotSingle == False:
                    fh_.write('ToKeepRow|1')
                    list_readsToKeep_1.append(contigName)


                fh_.write('\n')


    print("######### FOUND IN /2 only")
    with open(fn_only2, 'w+') as fh_:
        for contigName in dict_readsToIS_2:
            if contigName not in dict_readsToIS_1:
                fh_.write(contigName + '\t')
                for ISname in dict_readsToIS_2[contigName]:
                    for anInfo in dict_readsToIS_2[contigName][ISname]:
                        (alignLen, readLen, choppedLen) = getTheISAndInfo(anInfo)

                        isKeep="Rm"

                        if choppedLen >= minChoppedLen:
                            isKeep = "Keep"
                        else:
                            isAnyRm_2 = True


                        fh_.write(ISname + ":" + str(alignLen) + ":" + str(readLen) +  ":" + str(choppedLen) + ":" + isKeep +  ';')

                        arr_ISin2.append(ISname)

                fh_.write('\t')


                isAnyNotSingle = False
                if len(dict_readsToIS_2[contigName]) > 1:
                    fh_.write('Multiple')
                    isAnyNotSingle = True
                else:
                    fh_.write("Single")

                fh_.write('\t')
                if isKeep == "Keep" and isAnyNotSingle == False:
                    fh_.write('ToKeepRow|2')
                    list_readsToKeep_2.append(contigName)

                fh_.write('\n')

    return (list_readsToKeep_1, list_readsToKeep_2)


def loadTheReadNamesAndLen(fn_flanks, fn_out, list_complAlignedOther):

    fh_out = open(fn_out, 'w+')


    dict_readLens = dict() # dict_{readName} => len


    with open(fn_flanks, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            arr = record.id.split("__")

            if (arr[0] not in list_complAlignedOther):
                fh_out.write(">" + record.id + "\n")
                fh_out.write(str(record.seq) + "\n")
            else:
                print (arr[0])
            # print(arr[0])
            #if record.id not in list_complAlignedOther:
                # Print to file

            #    dict_readLens[record.id] = len(record.seq)

    # return (dict_readLens)
        

################################# MAIN
def main():
    parser = argparse.ArgumentParser(description='Load both reads files, for reads and their lengths, then check in blast results to see if either flank aligns 100%, then remove that from the flanks_1|2 files.')

    parser.add_argument('--blastRes_1', nargs=1, required=True, help="Blast results of reads to IS.")
    parser.add_argument('--blastRes_2', nargs=1, required=True, help="Blast results of reads to IS.")

    parser.add_argument('--minChoppedLen', nargs=1, default=[18], type=int, help="Minimun chopped sequence length (excluding overlapping sequence) to keep that sequence for further analysis. Default=18.")
    parser.add_argument("--buffer", nargs=1, default=[0], type=int, help="Buffer sequence length of the alignment of SKESA to complete genome to include, to see if the IS lies nearby. Default=0.")

    parser.add_argument('--flanks_1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--flanks_2', nargs=1, required=True, help="Flanks_2 file.")

    parser.add_argument('--out_flanks_1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--out_flanks_2', nargs=1, required=True, help="Flanks_2 file.")

    parser.add_argument('--flanks1AndOnly2', nargs=1, help="Output from summarizeOnlyToContig.py", default=[''])
    parser.add_argument('--flanks2AndOnly1', nargs=1, help="Output from summarizeOnlyToContig.py", default=[''])

    args = parser.parse_args()

    rmTheFlanks(args.blastRes_1[0], args.blastRes_2[0], args.minChoppedLen[0], args.flanks_1[0], args.flanks_2[0], args.out_flanks_1[0], args.out_flanks_2[0], args.flanks1AndOnly2[0], args.flanks2AndOnly1[0]) # args.flanks_1[0], args.flanks_2[0], args.out_flanks_1[0], args.out_flanks_2[0])

if __name__ == '__main__':
    main()
