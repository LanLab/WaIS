#!/usr/bin/python3 

import argparse
import re
import sys 
import os

def genPresenceAbsenceTable(fn_IStoRef_illuminaCalc, fn_IStoRef_blast):

	list_ISposFromContigs = load_IStoRef_illuminaCalc(fn_IStoRef_illuminaCalc)


	dict_refPresAbs = determinePresAbs(fn_IStoRef_blast, list_ISposFromContigs)
	
	# Print the output 
	sys.stdout.write(os.path.basename(fn_IStoRef_blast))
	for (start, end) in sorted(dict_refPresAbs): 
		sys.stdout.write ("\t" + str(start) + ":" + str(end)); 
	print(); 

	sys.stdout.write(os.path.basename(fn_IStoRef_illuminaCalc))
	for key in sorted(dict_refPresAbs): 
		sys.stdout.write ("\t" + str(dict_refPresAbs[key])); 
	print() 


def determinePresAbs(fn_IStoRef_blast, list_ISposFromContigs): 
	dict_refPos = dict() # dict_{(start, end)} => isPresent

	Col_refStart = 8
	Col_refEnd = 9
	with open (fn_IStoRef_blast, 'r') as fh: 
		for line in fh:  
			if re.match("^\#", line): 
				continue 

			arr = line.split("\t") 

			refStart = int(arr[Col_refStart])
			refEnd = int (arr[Col_refEnd])

			isPresent = checkIfPresent(refStart, refEnd, list_ISposFromContigs)
			# print (arr[Col_refStart] + '\t' + arr[Col_refEnd] + "\t" + str(isPresent))

			dict_refPos[(refStart, refEnd)] = isPresent

	return dict_refPos

def checkIfPresent(refStart, refEnd, list_ISposFromContigs):

	for (start_calcIS, end_calcIS) in list_ISposFromContigs: 

		if (refStart >= start_calcIS and refStart <= end_calcIS) or (refEnd >= start_calcIS and refEnd <= end_calcIS) or (start_calcIS >= refStart and start_calcIS <= refEnd) or (end_calcIS >= refStart and end_calcIS <= refEnd): 
			# print ("Found val: " + str(start_calcIS) + ":" + str(end_calcIS))
			return 1 

	return 0 
	

def load_IStoRef_illuminaCalc(fn_): 
	Col_start = 3 
	Col_end = 4 

	list_positions = [] # [(startPos, endPos), (startPos, endPos) ...]
	with open(fn_, 'r') as fh: 
		for line in fh: 
			if re.match("^\#", line): 
				continue

			line = line.strip() 
			arr = re.split("\t", line)

			# print (arr[Col_start] + '\t' + arr[Col_end])

			startPos = int(arr[Col_start]) 
			endPos = int(arr[Col_end])

			list_positions.append((startPos, endPos))

	return list_positions


#################################### MAIN 
def main(): 
	parser = argparse.ArgumentParser(description='Generate a presence absence table of IS w.r.t to a reference genome.')

	parser.add_argument('--IStoRef_illuminaCalc', nargs=1, required=True, help='IS positions in reference determined from contigs.')
	parser.add_argument('--IStoRef_blast', nargs=1, required=True, help='IS positions in reference determined from blast (in .gff format).')

	args = parser.parse_args()

	genPresenceAbsenceTable(args.IStoRef_illuminaCalc[0], args.IStoRef_blast[0])

if __name__ == '__main__':
	main() 
