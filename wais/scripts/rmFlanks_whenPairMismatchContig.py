#!/usr/bin/python3

import argparse
from Bio import SeqIO

#################################################### TOP_LVL
def filterFlanks(fn_only2Flanks1, fn_only1Flanks2, fn_flanks1, fn_flanks2, fn_out_flanks1, fn_out_flanks2):

    list_readIds_1 = loadFlankIdsToRm(fn_only2Flanks1)
    filterTheFlanks(fn_flanks1, list_readIds_1, fn_out_flanks1)


    list_readIds_2 = loadFlankIdsToRm(fn_only1Flanks2)
    filterTheFlanks(fn_flanks2, list_readIds_2, fn_out_flanks2)



#################################################### AUX
def filterTheFlanks(fn_flanks, list_readIdsToRm, fn_out_flanks):

    fh_out = open(fn_out_flanks, 'w+')

    with open (fn_flanks, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            if record.id not in list_readIdsToRm:
                fh_out.write(">" + record.id + "\n")
                fh_out.write(str(record.seq) + "\n")
    fh_out.close()


def loadFlankIdsToRm(fn_only2Flanks1):

    Col_readId = 0
    Col_isSame = 2
    list_readIds = []

    with open(fn_only2Flanks1, 'r') as fh:
        for line in fh:

            line = line.strip()
            arr = line.split("\t")

            if (len(arr) - 1 > Col_isSame) and arr[Col_isSame] == "False":
                # print (arr[Col_readId] + "\t" + arr[Col_isSame])
                list_readIds.append(arr[Col_readId])

    return list_readIds


#################################################### MAIN
def main():
    parser = argparse.ArgumentParser(description='Remove flanks whole complete pairs (i.e only files) ')

    parser.add_argument('--only1AndFlanks2', nargs=1, required=True, help="Do the pairs align to same contig.")
    parser.add_argument('--only2AndFlanks1', nargs=1, required=True, help="Do the pairs align to same contig.")

    parser.add_argument('--flanks1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--flanks2', nargs=1, required=True, help="Flanks_2 file.")

    parser.add_argument('--out_flanks1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--out_flanks2', nargs=1, required=True, help="Flanks_2 file.")


    args = parser.parse_args()

    filterFlanks(args.only2AndFlanks1[0], args.only1AndFlanks2[0], args.flanks1[0], args.flanks2[0], args.out_flanks1[0], args.out_flanks2[0])



if __name__ == '__main__':
    main()
