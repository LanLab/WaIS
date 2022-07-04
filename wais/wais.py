#!/usr/bin/python3 

import argparse
import sys
import os 

#################################### TOP_LVL 



#################################### AUX



#################################### MAIN 
def main(): 
	parser = argparse.ArgumentParser(description='Determine where the insertion sequences (IS) are in a reference genome, or an assembly (generated using short-reads) - using short-read sequences.')

	parser.add_argument('--outputDir', required=True, nargs=1)
	parser.add_argument('--runSpades', action='store_true', help='')
	parser.add_argument('--assembly', nargs=1, help='')
	parser.add_argument('--ISseqs', required=True, nargs=1, help='Fasta file containing the IS sequences to find.') 
	parser.add_argument('--reads_1', required=True, nargs=1, help="Illumina reads forward file.")
	parser.add_argument('--reads_2', required=True, nargs=1, help="Illumina reads reverse file.")



	## Spades options
	# --path_to_spades 
	# --options_to_spades

	## WaIS thresholds


	args = parser.parse_args()
	checkSpades(args.runSpades, args.assembly)
	checkIfFileExists(args.ISseqs[0])
	checkIfFileExists(args.reads_1[0])
	checkIfFileExists(args.reads_2[0])


def checkSpades(runSpades, assembly): 

	if ((runSpades == None or runSpades == False) and assembly == None) or (assembly and len(assembly) != 1): 
		sys.exit("\nError: atleast one of the following options is required:\n --runSpades\n --assembly <contigs.fasta>\n\n")

	if assembly: 
		checkIfFileExists(assembly[0])
		
def checkIfFileExists(filename): 		
	if not os.path.exists(filename): 
		sys.exit('\nError: the specified assembly file does not exist ' + filename + '\n')

if __name__=='__main__':
	main() 