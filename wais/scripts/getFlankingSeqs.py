import argparse
import re

from Bio import SeqIO


################################## TOP_LVL
def getTheFlanks(fn_blastRes_readsToIS, fn_reads, overlap, minChoppedLen):

    dict_parAlignReads = loadBlastRes(fn_blastRes_readsToIS)

    # dict_parAlignReads = dict()
    loadUnalignedSeqs(dict_parAlignReads, fn_reads, minChoppedLen)
    """
    for readName in dict_parAlignReads:
        for ISname in dict_parAlignReads[readName]:
            print (readName + "\t" + ISname + "\t" + str(dict_parAlignReads[readName][ISname]))
    pass
    """

################################## AUX
def loadUnalignedSeqs(dict_parAlignReads, fn_reads, minChoppedLen):
    sep = "__"

    with open(fn_reads, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):

            if record.id in dict_parAlignReads:

                for ISname in dict_parAlignReads[record.id]:
                    for (alignStart, alignEnd, readLen, unalignedSeq_pre, unalignedSeq_post, alignDir) in dict_parAlignReads[record.id][ISname]:

                    # (alignStart, alignEnd, readLen, unalignedSeq_pre, unalignedSeq_post) = dict_parAlignReads[record.id][ISname]

                        preStart = 0
                        preEnd = alignStart - 1

                        postStart = alignEnd -1
                        postEnd = readLen

                        if (preEnd - preStart) >= minChoppedLen:


                            # print (record.id + "\tPre\t" + ISname + "\t" + str(preStart) + "\t" + str(preEnd) + "\t" + str(postStart) + "\t" + str(postEnd))
                            print(">" + record.id + sep + "Pre" + sep + ISname + sep + str(preStart) + sep + str(preEnd) + sep + str(postStart) + sep + str(postEnd) + sep + alignDir)

                            print(record.seq[preStart:preEnd])
                        if (postEnd - postStart) >= minChoppedLen:

                            print (">" + record.id + sep +  "Post" + sep + ISname + sep + str(preStart) + sep + str(preEnd) + sep + str(postStart) + sep + str(postEnd) + sep + alignDir)

                            print(record.seq[postStart:postEnd])

                    # ExtractSeq. (also get the thresholded len??)





def loadBlastRes(fn_blastRes_readsToIS):

    dict_parAlignReads = dict() # Partially aligned reads # dict_{readNames} => {ISname} => (qStart, qEnd, qLen, qSeq_unalignedSegPre, qSeq_unalignedSeqPost, alignDir)

    Col_qSeqId = 0
    Col_sSeqId = 1 # Subject seq Id (ISname)
    Col_alignLen = 3
    Col_qStart = 6 # start of alignment in read
    Col_qEnd = 7 # end of alignment in read
    Col_sStart = 8 # start of alignment in IS
    Col_sEnd = 9 # end of alignment in IS
    Col_qLen = 12

    with open(fn_blastRes_readsToIS, 'r') as fh:
        for line in fh:

            if not re.match("^\#", line):

                line = line.strip()
                arr = line.split("\t")

                if int(arr[Col_alignLen]) < int(arr[Col_qLen]): # alignment length is less as read length;

                    # calcLenOfChopped(int(arr[Col_qStart]), int(arr[Col_qEnd]), int(arr[Col_qLen]));

                    alignDir = getAlignDir_readToIS(int(arr[Col_qStart]), int(arr[Col_qEnd]), int(arr[Col_sStart]), int(arr[Col_sEnd])) ;

                    if arr[Col_qSeqId] not in dict_parAlignReads:
                        dict_parAlignReads[arr[Col_qSeqId]] = dict()

                    if arr[Col_sSeqId] not in dict_parAlignReads[arr[Col_qSeqId]]:
                        dict_parAlignReads[arr[Col_qSeqId]][arr[Col_sSeqId]] = []

                    dict_parAlignReads[arr[Col_qSeqId]][arr[Col_sSeqId]].append((int(arr[Col_qStart]), int(arr[Col_qEnd]), int(arr[Col_qLen]), "", "", alignDir))

                    # print (arr[Col_qSeqId] + "\t" + arr[Col_alignLen] + "\t" +  arr[Col_qStart] + "\t" + arr[Col_qEnd] + "\t" + arr[Col_qLen])

    return (dict_parAlignReads)

def getAlignDir_readToIS(alignStart_read, alignEnd_read, alignStart_IS, alignEnd_IS):

    if (alignStart_IS <= alignEnd_IS):
        if (alignStart_read <= alignEnd_read):
            return '+'

        if (alignStart_read > alignEnd_read):
            return '-'

    if (alignStart_IS > alignEnd_IS):
        if (alignStart_read <= alignEnd_read):
            return '-'

        if (alignStart_read > alignEnd_read):
            return '+'

def calcLenOfChopped(readAlignStart, readAlignEnd, readLen):

    print (str(readAlignStart) + "\t" + str(readAlignEnd) + "\t" +  str(readLen))
    pass

################################## MAIN

def main():

    parser = argparse.ArgumentParser(description='Get reads which dont align to IS sequences (some IS-aligned-seq. can be added to the unaligned seqs - as skesa does contain alignments to IS)')

    parser.add_argument('--reads_to_IS', nargs=1, help="Blast results for reads to IS", required=True)
    parser.add_argument('--overlap', nargs=1, help="Number of base pairs of IS alignment to include", default=[0], type=int)
    parser.add_argument('--minChoppedLen', nargs=1, default=[18], type=int, help="Minimun chopped sequence length (excluding overlapping sequence) to keep that sequence for further analysis.")
    parser.add_argument('--reads', nargs=1, required=True, help="Reads filename; in fasta format.")

    args = parser.parse_args()

    getTheFlanks(args.reads_to_IS[0], args.reads[0], args.overlap[0], args.minChoppedLen[0])

if __name__ == '__main__':
    main()
