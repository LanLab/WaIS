#!/usr/bin/python3 

import argparse
import re
import sys 
import os

def genPresenceAbsenceTable(fn_IStoRef_illuminaCalc, fn_ISinRef_blast, isolateId):

	# 1. Load fn_IStoRef (the IS in the reference genome). 
	dict_ISinRef = load_ISinRef(fn_ISinRef_blast)

	# 2. Handle ISinRef_calculated 
	(dict_IStoRefCalc, dict_refIS_encountered) = handle_IStoRefCalc(fn_IStoRef_illuminaCalc)

	# 3. For each IStoRef, get ISinRefCalc (and note those found)  
	printPresenceAbsenceInRef(dict_ISinRef, dict_IStoRefCalc, dict_refIS_encountered, isolateId) 	

	"""
	for refId in dict_ISinRef: 
		for (refStart, refEnd) in dict_ISinRef[refId]: 
			print (refId + '\t' + str(refStart) + ':' + str(refEnd) + '\t' + str(dict_ISinRef[refId][(refStart, refEnd)]))
	""" 


"""
	list_ISposFromContigs = load_IStoRef_illuminaCalc(fn_IStoRef_illuminaCalc)

	print (list_ISposFromContigs)

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
"""
###################################### OUTPUT
def printPresenceAbsenceInRef(dict_ISinRef, dict_IStoRefCalc, dict_refIS_encountered, isolateId): 

	# print ('# IS in ref, also found in ' + isolateId)
	# Previously found in ref.
	isSecondOnwards = False 
	for refId in dict_ISinRef: 
		if isSecondOnwards == False:
			isSecondOnwards = True 
		else: 
			sys.stdout.write('\t') 
		sys.stdout.write(refId)

		for (refStart, refEnd) in sorted(dict_ISinRef[refId]): 
			for (ISid, orient) in dict_ISinRef[refId][(refStart, refEnd)]:

				refKey = str(refStart) + ':' + str(refEnd) + ':' + str(ISid) + ':' + str(orient)
				sys.stdout.write('\t' + refKey) 

	# sys.stdout.write('\n')

	# New, not previously found in ref.
	# isSecondOnwards = False
	for refId in dict_IStoRefCalc: 
		for idx, anIns in enumerate(dict_IStoRefCalc[refId]['Insertions']):
			if dict_IStoRefCalc[refId]['isFoundInRef'][idx] == False: 
				sys.stdout.write ('\t' + str(anIns)) 
			# sys.stdout.write (dict_IStoRefCalc[refId]['isFoundInRef'])

	sys.stdout.write('\n')

	############################### Printing for ref

	isSecondOnwards = False 
	for refId in dict_ISinRef: 
		if isSecondOnwards == False:
			isSecondOnwards = True 
		else: 
			sys.stdout.write('\t') 
		sys.stdout.write(refId)

		for (refStart, refEnd) in sorted(dict_ISinRef[refId]): 
			for (ISid, orient) in dict_ISinRef[refId][(refStart, refEnd)]:

				# refKey = str(refStart) + ':' + str(refEnd) + ':' + str(ISid) + ':' + str(orient)
				sys.stdout.write('\t' + '1') 

	# sys.stdout.write('\n')

	# New, not previously found in ref.
	# isSecondOnwards = False
	for refId in dict_IStoRefCalc: 
		for idx, anIns in enumerate(dict_IStoRefCalc[refId]['Insertions']):
			if dict_IStoRefCalc[refId]['isFoundInRef'][idx] == False: 
				sys.stdout.write ('\t' + '0') 
			# sys.stdout.write (dict_IStoRefCalc[refId]['isFoundInRef'])

	sys.stdout.write('\n')



	############################### Printing for isolate
				
	sys.stdout.write(isolateId)

	isSecondOnwards = False 
	for refId in dict_ISinRef: 
		if isSecondOnwards == False:
			isSecondOnwards = True 
		else: 
			sys.stdout.write('\t') 

		for (refStart, refEnd) in sorted(dict_ISinRef[refId]): 
			for (ISid, orient) in dict_ISinRef[refId][(refStart, refEnd)]:


				sys.stdout.write('\t')
				refKey = str(refStart) + ':' + str(refEnd) + ':' + str(ISid) + ':' + str(orient)
				
				if refId in dict_refIS_encountered and refKey in dict_refIS_encountered[refId]: 
					sys.stdout.write('1')
				else: 
					sys.stdout.write('0')

					

	# sys.stdout.write('\n')	

	
	# print ('# IS in ' + isolateId + ' not previously found in ref.')

	for refId in dict_IStoRefCalc: 		
		# print (len (dict_IStoRefCalc[refId]['isFoundInRef']))
		# print (len (dict_IStoRefCalc[refId]['Insertions']))
		for idx, anIns in enumerate(dict_IStoRefCalc[refId]['Insertions']):
			if dict_IStoRefCalc[refId]['isFoundInRef'][idx] == False: 
				sys.stdout.write ('\t' + '1') 
			# sys.stdout.write (dict_IStoRefCalc[refId]['isFoundInRef'])

	sys.stdout.write('\n')
	 

###################################### AUX
def handle_IStoRefCalc(fn_IStoRefCalc):

	dict_IStoRefCalc = dict() # dict_[refId] => {'Insertions': [(refStart, refEnd, ISid, orient), (), ... ], 'isFoundInRef': [T|F, T|F, T|F, ...]}

	dict_refIS_encountered = dict() # dict_[refId] => list_["refStart:refEnd:ISid:orient"]

	Col_refId = 0 
	Col_refStart = 3
	Col_refEnd = 4 
	Col_orient = 6
	Col_info = 8 

	with open(fn_IStoRefCalc, 'r') as fh: 
		for line in fh: 

			if not re.match(r'^\#', line): 
				line = line.strip()
				arr = line.split('\t') 
				
				(ISid, refIS) = extractISidFromInfo(arr[Col_info])
				# print (arr[Col_refId] + '\t' + arr[Col_refStart] + ':' + arr[Col_refEnd] + '\t' + arr[Col_orient] + '\t' + arr[Col_info] + '\t' + ISid) 

				if arr[Col_refId] not in dict_IStoRefCalc: 
					dict_IStoRefCalc[arr[Col_refId]] = dict()
					dict_IStoRefCalc[arr[Col_refId]]['Insertions'] = list() 
					dict_IStoRefCalc[arr[Col_refId]]['isFoundInRef'] = list() 

					dict_refIS_encountered[arr[Col_refId]] = list() 
					
				theInsertion = (int(arr[Col_refStart]), int(arr[Col_refEnd]), ISid, arr[Col_orient], refIS)
				if theInsertion not in dict_IStoRefCalc[arr[Col_refId]]: 
					dict_IStoRefCalc[arr[Col_refId]]['Insertions'].append(theInsertion)

					if refIS != '': 
						if ',' in refIS: 
							arr_refIS = refIS.split(',')
						else: 
							arr_refIS = [refIS]

						for anRefIS in arr_refIS: 
							if anRefIS not in dict_refIS_encountered[arr[Col_refId]]: 
							# 	dict_IStoRefCalc[arr[Col_refId]]['isFoundInRef'].append(True)
								dict_refIS_encountered[arr[Col_refId]].append(anRefIS)

							# else: 
						dict_IStoRefCalc[arr[Col_refId]]['isFoundInRef'].append(True)
					else: 
						dict_IStoRefCalc[arr[Col_refId]]['isFoundInRef'].append(False)

	return (dict_IStoRefCalc, dict_refIS_encountered)
				

def extractISidFromInfo(info): 

	arr = info.split(';')
	ISid = '' 
	refIS = '' 

	for val in arr: 
		if re.match('^Name\=', val): 
			ISid = re.sub('^Name\=', '', val)
			# print (val) 
		if re.match('^refIS\=', val): 
			refIS = re.sub('^refIS\=', '', val)
			# print (refIS)

		
	# print (refIS)
	if ISid == '': 
		sys.exit('Name= (i.e. the ISid) not found in the GFF file')

	return (ISid, refIS)  


def load_ISinRef(fn_): 

	dict_ISinRef = dict() # dict_[refId] => {(start, end)} => [(ISid, orient)]

	Col_refStart = 8 
	Col_refEnd = 9 
	Col_ISstart = 6 
	Col_ISend = 7 
	Col_ISid = 0 
	Col_refId = 1 


	with open(fn_, 'r') as fh: 
		for line in fh: 

			if not re.match(r'^\#', line): 
				line = line.strip() 

				arr = line.split('\t')

				
				# print (line)  
				# print (arr[Col_refId] + '\t' + arr[Col_ISid] + '\t' + arr[Col_refStart] + ":" + arr[Col_refEnd] + '\t' + arr[Col_ISstart] + ':' + arr[Col_ISend])

				addToDict_ISinRef(dict_ISinRef, arr[Col_refId], arr[Col_ISid], int(arr[Col_refStart]), int(arr[Col_refEnd]), int(arr[Col_ISstart]), int(arr[Col_ISend])) 

	
	return dict_ISinRef

def addToDict_ISinRef(dict_, refId, ISid, refStart, refEnd, ISstart, ISend): 

	if (refStart < refEnd and ISstart < ISend) or (refStart > refEnd and ISstart > ISend): 
		orient = '+'
	else: 
		orient = '-'

	if refStart > refEnd: 
		tmp = refStart 
		refStart = refEnd 
		refEnd = tmp 


	if refId not in dict_: 
		dict_[refId] = dict() 

	if (refStart, refEnd) not in dict_[refId]: 
		dict_[refId][(refStart, refEnd)] = list() 

	if (ISid, orient) not in dict_[refId][(refStart, refEnd)]: 
		dict_[refId][(refStart, refEnd)].append((ISid, orient))


	# print ('Calculated orientation: ' + orient + ' ' + str(refStart) + ':' + str(refEnd) + ' ' + str(ISstart) + ':' + str(ISend))
###################################### AUX - PREV 

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

	parser.add_argument('--IStoRef_mapped', nargs=1, required=True, help='IS positions in reference (in GFF3 format) determined from contigs using the short reads.')
	parser.add_argument('--IStoRef_blast', nargs=1, required=True, help='IS positions in reference (in blast --outfmt7) determined from blast.')
	# parser.add_argument('--isMerged', action='store_true', required=False, default=False, help='Is provided file (IStoRef_mapped) just at the merged level')
	# parser.add_argument('--isIgnoreOrient', action='store_true', required=False, default=False, help='Is provided file (IStoRef_mapped) just at the ignoreOrient level')
	# parser.add_argument('--isIgnoreIStype', action='store_true', required=False, default=False, help='Is provided file (IStoRef_mapped) just at the ignoreIStype level')

	parser.add_argument('--isolateId', nargs=1, required=True, help='Isolate id, for printing in the output file.')
	


	args = parser.parse_args()

	# if (int(args.isMerged) + int(args.isIgnoreOrient) + int(args.isIgnoreIStype)) != 1: 
	# 	sys.exit("Error: exactly one of --isMerged, --isIgnoreOrient or --isIgnoreIStype should be specified. (Just choose the specific level - thank you :).")

	# print (args.isMerged + '\t' + args.isIgnoreOrient, args.isIgnoreIStype)

	genPresenceAbsenceTable(args.IStoRef_mapped[0], args.IStoRef_blast[0], args.isolateId[0])

if __name__ == '__main__':
	main() 
