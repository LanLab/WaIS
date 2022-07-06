#!/usr/bin/python3

import argparse
import re

from Bio import SeqIO

########################################### TOP_LVL
def loadAndPrint(fn_only1, fn_only2, fn_reads1, fn_reads2, fn_out_reads1, fn_out_reads2):

    list_readIds_1 = loadOnlyFile(fn_only1)
    list_readIds_2 = loadOnlyFile(fn_only2)

    printReadsToBlast(fn_reads2, list_readIds_1, fn_out_reads2)
    printReadsToBlast(fn_reads1, list_readIds_2, fn_out_reads1)

########################################### AUX
def printReadsToBlast(fn_reads, list_readIds, fn_out):
    # sep = "__"

    fh_out = open(fn_out, 'w+');

    with open(fn_reads, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            # print (record.id)

            if (record.id in list_readIds):
                fh_out.write ('>' + record.id + '\n')
                fh_out.write (str(record.seq) + '\n')


def loadOnlyFile(fn_only1):
    Col_readId = 0
    Col_isSingle = 2
    Col_isToKeepRow = 3

    list_readIds = []

    with open(fn_only1, 'r') as fh:
        for line in fh:
            line = line.strip()

            arr = line.split("\t")

            #print (len(arr))
            if (Col_isToKeepRow == (len(arr) -1) and arr[Col_isSingle] == 'Single' and re.match('ToKeepRow', arr[Col_isToKeepRow])):
                # print(arr[Col_readId] + "\t" + arr[Col_isSingle] + "\t" + arr[Col_isToKeepRow])
                list_readIds.append(arr[Col_readId])

    return list_readIds

########################################### MAIN
def main():
    parser = argparse.ArgumentParser(description='Load only1, and only2 files, and print out the fasta files where ToKeepRow and Single.')

    parser.add_argument('--only1', nargs=1, required=True, help="IS in only1.txt")
    parser.add_argument('--only2', nargs=1, required=True, help="IS in only2.txt")

    parser.add_argument('--reads1', nargs=1, required=True, help="Reads file")
    parser.add_argument('--reads2', nargs=1, required=True, help="Reads file")

    parser.add_argument('--outfile_reads1', nargs=1, required=True, help="Fasta outfile reads1")
    parser.add_argument('--outfile_reads2', nargs=1, required=True, help="Fasta outfile reads2")

    """
    parser.add_argument('--minChoppedLen', nargs=1, default=[18], type=int, help="Minimun chopped sequence length (excluding overlapping sequence) to keep that sequence for further analysis.")
    parser.add_argument("--buffer", nargs=1, default=[0], type=int, help="Buffer sequence length of the alignment of SKESA to complete genome to include, to see if the IS lies nearby.")

    parser.add_argument('--flanks_1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--flanks_2', nargs=1, required=True, help="Flanks_2 file.")

    parser.add_argument('--out_flanks_1', nargs=1, required=True, help="Flanks_1 file.")
    parser.add_argument('--out_flanks_2', nargs=1, required=True, help="Flanks_2 file.")
    """

    args = parser.parse_args()

    loadAndPrint(args.only1[0], args.only2[0], args.reads1[0], args.reads2[0], args.outfile_reads1[0], args.outfile_reads2[0])


if __name__ == '__main__':
    main()
