#!/usr/bin/python3

from calcInRef_posISorient_v2 import isAnyOverlap 
import calcInContig_posISOrient as calcInContig
import argparse
import sys
import re

SO_ = 'insertion_serquence'

#################################### TOP-LVL 
def convertToGff(fn_blastRes, fn_out, isMerged, ignoreOrient, ignoreIStype, separator, fn_estimates, fn_estimates_singleRow): 

	print ('isMerged ' + str(isMerged) + ' ' + str(ignoreOrient) + ' ' + str(ignoreIStype))

	(dict_ISinRef, refsEncountered) = loadISinRef(fn_blastRes) 

	"""
	for key in dict_ISinRef: 
		print (key) 
		print (dict_ISinRef[key])
	"""

	if isMerged == True: 
		dict_merged = checkAllTheISinRef(dict_ISinRef, ignoreOrient, ignoreIStype, separator)
		printGff3File(fn_out, refsEncountered, dict_merged) 
		dict_counts = calculateEstimates(dict_merged)

		# for fn_estimates in list_fn_estimates: 
		# 	print ('Writing estimates to ' + fn_estimates)
		printEstimates(fn_estimates, dict_counts) 
		printEstimates_singleRow(fn_estimates_singleRow, dict_counts)

#################################### OUTPUT
def printEstimates_singleRow(fn_estimates_singleRow, dict_counts):
	lines = None 

	with open(fn_estimates_singleRow, 'r') as fh:
		lines = fh.read()

	with open (fn_estimates_singleRow, 'w+') as fh: 
		for idx, line in enumerate(lines.split('\n')): 
			if idx == 0: 
				# Header line
				line = line + '\t' + '#Reference'

				for refId in sorted(dict_counts):
					line = line + '\t' + refId
					for ISid in sorted(dict_counts[refId]):
						line = line + '\t' + ISid
					line = line + '\t' + 'Total'
					 
				fh.write(line + '\n')

			if idx == 1: 
				line = line + '\t'

				total = 0
				for refId in sorted(dict_counts):
					line = line + '\t'
					for ISid in sorted(dict_counts[refId]):
						line = line + '\t' + str(dict_counts[refId][ISid]) 
						total = total + dict_counts[refId][ISid]
					line = line + '\t' + str(total)


				# line = line + '\t' + str(total)
				fh.write(line + '\n')
				pass 	

			# print (idx); 
			# print (line);

		# print(lines)
	 

def printEstimates(fn_estimates, dict_counts):

	with open(fn_estimates, 'a+') as fh: 

		fh.write('#################### Reference' + '\n')

		for refId in dict_counts:
			fh.write(refId) 
			for ISid in sorted(dict_counts[refId]):
				fh.write('\t' + ISid) 
			fh.write('\t' + 'Total' + '\n')

		total = 0 
		for refId in dict_counts:
			
			for ISid in sorted(dict_counts[refId]):
				fh.write('\t' + str(dict_counts[refId][ISid])) 
				total = total + dict_counts[refId][ISid]
			fh.write('\t' + str(total) + '\n')
	

		


def printGff3File(fn_out, refsEncountered, dict_merged):

	with open(fn_out, 'w+') as fh: 

		fh.write('##gff-version 3.1.26' + '\n')
		
		
		for refId in dict_merged: 
			for (encountered_refId, encountered_len) in refsEncountered: 
				if encountered_refId == refId: 
					fh.write('##sequence-region ' + refId + ' ' + '1' + ' ' + str(encountered_len) + '\n')

		for refId in dict_merged: 
			for (refStart, refEnd, refIS, refOrient) in dict_merged[refId]: 
				calcInContig.printGff3Line(fh, refId, 'BLAST+', SO_, refStart, refEnd, '.', refOrient, '.', 'Name=' + refIS)

		

#################################### AUX 
def calculateEstimates(dict_merged): 

	dict_counts = dict() # {refId} => {ISid} => count
	for refId in dict_merged: 

		if refId not in dict_counts: 
			dict_counts[refId] = dict() 
		
		for (refStart, refEnd, ISid, orient) in dict_merged[refId]: 

			if ISid not in dict_counts[refId]: 
				dict_counts[refId][ISid] = 0

			dict_counts[refId][ISid] = dict_counts[refId][ISid] + 1 
			
	return dict_counts 

def checkAllTheISinRef(dict_ISinRef, ignoreOrient, ignoreIStype, separator):

	dict_merged = dict() # {refId} => {(start, end, ISid, orient)} 
	for refId in dict_ISinRef: 
		for (refStart, refEnd, ISid, orient) in dict_ISinRef[refId]:
			checkIfNeedToBeMerged(dict_merged, refId, refStart, refEnd, ISid, orient, ignoreOrient, ignoreIStype, separator) 

	return dict_merged		

def checkIfNeedToBeMerged(dict_merged, refId, refStart, refEnd, refIS, refOrient, ignoreOrient, ignoreIStype, separator): 

	if refId not in dict_merged: 
		# add to dict_merged 
		dict_merged[refId] = list() 
		dict_merged[refId].append((refStart, refEnd, refIS, refOrient))

		return  

	isFound = False; idxToDel = -1; refIdOfInterest = ''
	newStart = -1 
	newEnd = -1
	newISid = '' 
	newOrient = ''

	
	for mergedRefId in dict_merged: 
		# print (mergedRefId)

		for idx, (mergedStart, mergedEnd, mergedIS, mergedOrient) in enumerate(dict_merged[mergedRefId]):

			if isFound == False and  refIS == mergedIS and refOrient == mergedOrient and isAnyOverlap(refStart, refEnd, mergedStart, mergedEnd, 0):
				
				newStart = refStart if refStart < mergedStart else mergedStart 
				newEnd = refEnd if refEnd > mergedEnd else mergedEnd
				newISid = refIS 
				newOrient = refOrient 

				idxToDel =  idx
				refIdOfInterest = mergedRefId
				isFound = True 

			elif isFound == False and ignoreIStype == True and isAnyOverlap(refStart, refEnd, mergedStart, mergedEnd, 0):

				newStart = refStart if refStart < mergedStart else mergedStart 
				newEnd = refEnd if refEnd > mergedEnd else mergedEnd
				
				newOrient = mergedOrient.split(separator)
				if refOrient not in newOrient: 
					newOrient.append(refOrient)
				
				newOrient = separator.join(sorted(newOrient)) 



				newISid = mergedIS.split(separator) 
				if refIS not in newISid: 
					newISid.append(refIS)

				newISid = separator.join(sorted(newISid))	
				

				idxToDel =  idx
				refIdOfInterest = mergedRefId
				isFound = True 	

			elif isFound == False and ignoreOrient == True and refIS == mergedIS and isAnyOverlap(refStart, refEnd, mergedStart, mergedEnd, 0): 

				# print ('This situation is happening ' + refIS + ' ' + mergedIS)

				newStart = refStart if refStart < mergedStart else mergedStart 
				newEnd = refEnd if refEnd > mergedEnd else mergedEnd
				
				newOrient = mergedOrient.split(separator)
				if refOrient not in newOrient: 
					newOrient.append(refOrient)
				
				newOrient = separator.join(newOrient) 
				newISid = refIS
				

				idxToDel =  idx
				refIdOfInterest = mergedRefId
				isFound = True 




	if isFound == True:
		del dict_merged[refIdOfInterest][idxToDel]
		dict_merged[refIdOfInterest].append((newStart, newEnd, newISid, newOrient))
	else: 
		dict_merged[refId].append((refStart, refEnd, refIS, refOrient))




def loadISinRef(fn_blastRes): 
	# list_ISinRef_loaded = list() 


	dict_ISinRef = dict() # dict_{refId} => {(refStart, refEnd)} => [(ISid, orient), ...]

	# fh_out = open(fn_out, 'w+')

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

				if arr[Col_sId] not in dict_ISinRef: 
					dict_ISinRef[arr[Col_sId]] = list() 
				if (refStart, refEnd, arr[Col_qId], orient) not in dict_ISinRef[arr[Col_sId]]: 
					dict_ISinRef[arr[Col_sId]].append((refStart, refEnd, arr[Col_qId], orient))
			


			   #  attributes = 'Name=' + arr[Col_qId]
				# calcInContig.printGff3Line(fh_out, arr[Col_sId], 'BLAST+', SO_, refStart, refEnd, '.', orient, '.', attributes)

				# list_ISinRef_loaded.append((arr[Col_sId], 'BLAST+', SO_, refStart, refEnd, '.', orient, '.', attributes))

				# print (line) 

	# fh_out.close()  


	###### 
	"""
	with open(fn_out, 'r') as original: 
		data = original.read()
	
	with open(fn_out, 'w') as modified: 
		modified.write('##gff-version 3.1.26' + '\n')
		for (contigId, contigLen) in refsEncountered: 
			modified.write('##sequence-region ' + contigId + ' ' + '1' + ' ' + str(contigLen) + '\n')
		
		modified.write(data)
	"""
	return (dict_ISinRef, refsEncountered)



#################################### MAIN 
def main(): 
	parser = argparse.ArgumentParser(description='Convert IS-to-reference blast results to gff. Exact blast alignment numbers are appended to the supplied estimates.txt file.')

	parser.add_argument('--blastRes', nargs=1, required=True, help='Blast results of contigs to reference.')
	parser.add_argument('--out', nargs=1, required=True, help='Output filename.')

	parser.add_argument('--isMerged', nargs=1, type=bool, default=[True], help='Make sure outputs are unique for a set of positions, IStype and orientation. Note:TODO: this option set to \'False\' is not programmed (there will be no output GFF3 file).')
	parser.add_argument('--ignoreOrient', nargs=1, type=bool, default=[False], help='Outputs are printed such that orientation is ignored.')
	parser.add_argument('--ignoreIStype', nargs=1, type=bool, default=[False], help='Outputs are printed such that IStype (and orientation) is ignored.')

	parser.add_argument('--separator', nargs=1, type=bool, default=[':'], help='Separator to apply between two different orientations or IStypes.')

	parser.add_argument('--fn_estimates', nargs='+', required=True, help='An estimates.txt file to which the reference IS insertions are appended.')
	parser.add_argument('--fn_estimates_singleRow', nargs='+', required=True, help='An estimates.txt file (in single row) to which the reference IS insertions are appended.')
	
	args = parser.parse_args()

	# print (args.fn_estimates)

	convertToGff(args.blastRes[0], args.out[0], bool(args.isMerged[0]), bool(args.ignoreOrient[0]), bool(args.ignoreIStype[0]), args.separator[0], args.fn_estimates[0], args.fn_estimates_singleRow[0])

if __name__ == '__main__':
	main()