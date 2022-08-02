#!/usr/bin/python3

import calcInContig_posISOrient as calcInContig
import re
import argparse

SO_ = 'contig'

#################################### TOP-LVL 
def convertToGff(fn_blastRes, fn_out): 

    fh_out = open(fn_out, 'w+')

    refsEncountered = [] 

    Col_qId = 0 
    Col_sId = 1 
    Col_sStart = 8
    Col_sEnd = 9 
    
    Col_qStart = 6
    Col_qEnd = 7

    Col_qlen = 12
    Col_slen = 13
    with open (fn_blastRes, 'r') as fh: 
        for line in fh: 
            
            if not re.match('^#', line): 
                line = line.strip() 
                arr = line.split('\t') 

                calcInContig.addTo_contigsEncountered(refsEncountered, (arr[Col_sId], arr[Col_slen]))

                refStart = int(arr[Col_sStart])
                refEnd = int(arr[Col_sEnd])

                contigStart = int(arr[Col_qStart])
                contigEnd = int(arr[Col_qEnd])
                
                orient = '.'
                if (refStart > refEnd and contigStart > contigEnd) or (refStart < refEnd and contigStart < contigEnd): 
                    orient = '+'
                else: 
                    orient = '-'


                if refStart > refEnd: 
                    tmp = refStart
                    refStart = refEnd
                    refEnd = tmp 

                attributes = 'Name=' + arr[Col_qId]
                calcInContig.printGff3Line(fh_out, arr[Col_sId], 'BLAST+', SO_, refStart, refEnd, '.', orient, '.', attributes)
                print (line) 

    fh_out.close()  


    with open(fn_out, 'r') as original: 
        data = original.read()
    
    with open(fn_out, 'w') as modified: 
        modified.write('##gff-version 3.1.26' + '\n')
        for (contigId, contigLen) in refsEncountered: 
            modified.write('##sequence-region ' + contigId + ' ' + '1' + ' ' + str(contigLen) + '\n')
        
        modified.write(data)

#################################### AUX 

#################################### MAIN 
def main(): 
    parser = argparse.ArgumentParser(description='Convert contigs-to-reference blast results to gff.')



    parser.add_argument('--blastRes', nargs=1, required=True, help='Blast results of contigs to reference.')
    parser.add_argument('--out', nargs=1, required=True, help='Output filename.')


    args = parser.parse_args()

    convertToGff(args.blastRes[0], args.out[0])

if __name__ == '__main__':
    main()