#!/usr/bin/python3 

import argparse 
from BCBio import GFF 
from Bio import SeqIO 
import sys
import re

UPSTREAM = 'Upstream'
DOWNSTREAM = 'Downstream'
DIRECT = 'Interruption'

################################################## TOP_LVL 
def insertionsWrtRefAnnotations(fn_reference, fn_refAnnotation, fn_ISannotation, allowedFeatTypes, th_surroundDist, fn_outfile_all, fn_outfile_onlyInterruped):

	(dict_refFeatures, dict_refSeq) = load_refAnnot(fn_reference, fn_refAnnotation, allowedFeatTypes) 

	"""
	for recordId in dict_refFeatures: 
		for (start, end, orient) in dict_refFeatures[recordId]: 
			print (recordId + ' ' + str(start) + ' ' + str(end) + ' ' + orient)
	"""
	
	dict_ISannots = load_ISannot(fn_ISannotation, dict_refSeq)
	
	"""
	# print (dict_ISannots.keys())
	for recordId in dict_ISannots:  
		for (start, end, orient) in dict_ISannots[recordId]: 
			print (recordId + ' ' + str(start) + ':' + str(end) + ':' + orient + ' ' + str(dict_ISannots[recordId][(start, end, orient)]))
	"""
	
	
	dict_refAnnots_interrupted = getInterruptedAnnots(dict_refFeatures, dict_ISannots) # Annotations directly disrupted by IS (i.e. not the surrounding region).

	
	"""
	for refId in dict_refAnnots_interrupted: 
		for refInfo in dict_refAnnots_interrupted[refId]: 

			print (refId + '\t' + str(refInfo) + '\t' + str(dict_refAnnots_interrupted[refId][refInfo]))
	"""
	
	 
	if th_surroundDist > 0: 
		# Calc effective distance to check
		dict_distToCheckForInsert = calcEffectiveDistToCheck(dict_refFeatures, th_surroundDist)

		dict_insertsInDist = getInterruptionsInDist(dict_distToCheckForInsert, dict_ISannots, dict_refAnnots_interrupted)

		
		#for refId in dict_distToCheckForInsert: 
		#	for (refStart, refEnd, refOrient) in dict_distToCheckForInsert[refId]:
		#		print (str(refStart) + ':' + str(refEnd) + ':' + str(refOrient) + '\t' + str(dict_distToCheckForInsert[refId][(refStart, refEnd, refOrient)])) 

	
	
		printFinalOutput(dict_refFeatures, dict_ISannots, dict_refAnnots_interrupted, th_surroundDist, dict_distToCheckForInsert, dict_insertsInDist, fn_outfile_all)

		if fn_outfile_onlyInterruped: 
			printFinalOutput_onlyInterrupted(fn_outfile_all, True, fn_outfile_onlyInterruped)


	else: 
		printFinalOutput(dict_refFeatures, dict_ISannots, dict_refAnnots_interrupted, th_surroundDist, None, None, fn_outfile_all)
 
		if fn_outfile_onlyInterruped: 
			printFinalOutput_onlyInterrupted(fn_outfile_all, False, fn_outfile_onlyInterruped)
	


	# dict_distToCheckForInsert

	"""
	for key in dict_refAnnots_interrupted: 
		print (str(key) + '\t' + str(dict_refAnnots_interrupted[key])) 


	for key in dict_refFeatures: 
		print (key) 
		for feature in dict_refFeatures[key]: 
			print (feature)
	"""
	# print (dict_refFeatures)
 
################################################## AUX - OUTPUT 
def printFinalOutput_onlyInterrupted(fn_outfile_all, isLeftRight, fn_outfile_onlyInterruped): 

	fh_out = open(fn_outfile_onlyInterruped, 'w+')
	fh_out.write('Ref_id' + '\t' + 'Ref_loc' + '\t' + 'Ref_annot_details' + '\t' + 'Interruption_type' + '\t' + 'Interrupted_by' + '\t' + 'IsISinsertionNew' + '\t' + 'InsertionPos_or_DistanceFromAnnot' + '\t' + 'EffectiveDistanceChecked')
	fh_out.write('\n')


	isHeader = True 

	Col_refId = 0 
	Col_refLoc = 1 
	Col_refAnnotDetails = 2 
	Col_isInterrupted = 3
	Col_interruptedBy = 4
	NumElems_interruptedBy = 5
	
	Col_refLoc_start = 0
	Col_refLoc_end = 1
	Col_refLoc_orient = 2
	

	Col_surroundDist = 5
	Col_isLeftInterup = 6 
	Col_left_interuptBy = 7 
	Col_isRightInterup = 8 
	Col_right_interupBy = 9

	# Each inserting IS is presented one per line.

	with open (fn_outfile_all, 'r') as fh: 
		for line in fh: 
			
			if isHeader == True: 
				isHeader = False
				continue 
				
			line = line.strip() 

			arr = line.split('\t')

			arr_refLoc = arr[Col_refLoc].split(':')

			# print ('Interrupted by: ' + arr[Col_interruptedBy])

			if arr[Col_isInterrupted] == 'True': 
				printAnInsertionLine(DIRECT, fh_out, arr[Col_refId], arr[Col_refLoc], arr[Col_refAnnotDetails], arr[Col_interruptedBy], '')


			if isLeftRight == True and arr[Col_isLeftInterup] == 'True':
				if arr_refLoc[Col_refLoc_orient] == '+': 
					# print ('It is upstream ' + str(arr[Col_refLoc]) + ' ' + arr_refLoc[Col_refLoc_orient]) 
					(leftDistChecked, rightDistChecked) = arr[Col_surroundDist].split(',')
					printAnInsertionLine(UPSTREAM, fh_out, arr[Col_refId], arr[Col_refLoc], arr[Col_refAnnotDetails], arr[Col_left_interuptBy], leftDistChecked)

					# printUpstreamDownstreamLines(UPSTREAM, fh, )

				elif arr_refLoc[Col_refLoc_orient] == '-':
					# print ('It is downstream ' + str(arr[Col_refLoc]) + ' ' + arr_refLoc[Col_refLoc_orient])  
					(leftDistChecked, rightDistChecked) = arr[Col_surroundDist].split(',')
					printAnInsertionLine(DOWNSTREAM, fh_out, arr[Col_refId], arr[Col_refLoc], arr[Col_refAnnotDetails], arr[Col_left_interuptBy], leftDistChecked)
					# printUpstreamDownstreamLines(DOWNSTREAM, fh, )

			if isLeftRight == True and  arr[Col_isRightInterup] == 'True':
				if arr_refLoc[Col_refLoc_orient] == '+': 
					# print ('It is downstream ' + str(arr[Col_refLoc]) + ' ' + arr_refLoc[Col_refLoc_orient]) 
					(leftDistChecked, rightDistChecked) = arr[Col_surroundDist].split(',')
					printAnInsertionLine(DOWNSTREAM, fh_out, arr[Col_refId], arr[Col_refLoc], arr[Col_refAnnotDetails], arr[Col_right_interupBy], rightDistChecked)
					# printUpstreamDownstreamLines(DOWNSTREAM, fh, )

				elif arr_refLoc[Col_refLoc_orient] == '-':
					# print ('It is upstream ' + str(arr[Col_refLoc]) + ' ' + arr_refLoc[Col_refLoc_orient])  
					(leftDistChecked, rightDistChecked) = arr[Col_surroundDist].split(',')
					printAnInsertionLine(UPSTREAM, fh_out, arr[Col_refId], arr[Col_refLoc], arr[Col_refAnnotDetails], arr[Col_right_interupBy], rightDistChecked)
					# printUpstreamDownstreamLines(UPSTREAM, fh)
			
				#
				# print (arr_interruptedByIS)
	 

def printAnInsertionLine(insertion_location, fh_out, refId, refLoc, refAnnotDetails, interruptingIS, distChecked):

	arr_interruptedByIS = interruptingIS.split(';')

	#  print (arr_interruptedByIS)

	for anInterruptingIS in arr_interruptedByIS: 
		# print ('An interrupting IS ' + anInterruptingIS)
		if ':' in anInterruptingIS: # len(arr_anInterup) == NumElems_interruptedBy:
			# print (arr_anInterup)

			arr_interupIS_elems = anInterruptingIS.split(':')
			insertPos_or_dist = arr_interupIS_elems.pop()

			insertPos_or_dist = re.sub('insertPos=', '', insertPos_or_dist)
			insertPos_or_dist = re.sub('leftDist=', '', insertPos_or_dist)
			insertPos_or_dist = re.sub('rightDist=', '', insertPos_or_dist)
			# insertPos_or_dist = re.sub('insertPos=', '', insertPos_or_dist)

			
			isNewInsertion = arr_interupIS_elems.pop()
			isNewInsertion = re.sub('isNew=', '', isNewInsertion)


			fh_out.write(refId + '\t' + refLoc + '\t' + refAnnotDetails + '\t' + insertion_location + '\t' + ":".join(arr_interupIS_elems) + '\t' + isNewInsertion + '\t' + insertPos_or_dist  + '\t' + distChecked)
			fh_out.write('\n')
	

def printFinalOutput(dict_refFeatures, dict_ISannots, dict_refAnnots_interrupted, th_surroundDist, dict_distToCheckForInsert, dict_insertsInDist, fn_outfile): 

	try: 
		fh_out = open(fn_outfile, 'w+')
	except: 
		fh_out.write("Error: could not open " + fn_outfile + '. Please check the path and try again.\n')
		sys.exit() 
	# Printing header
	fh_out.write ('Ref_id' + '\t' + 'Ref_loc' + '\t' + 'Ref_annot_details' + '\t' + 'isInterrupted' + '\t' + 'InterruptedBy') # Interruption position in reference #TODO

	if th_surroundDist > 0: 
		fh_out.write('\t' + 'Effective_surround_dist(left),(right)')

		fh_out.write('\t' + 'isLeftInterrupted')
		fh_out.write('\t' + 'Left_interruptedBy')
		
		fh_out.write('\t' + 'isRightInterrupted')
		fh_out.write('\t' + 'Right_interruptedBy') 


	fh_out.write ('\n')

	for refId in dict_refFeatures: 
		for (refStart, refEnd, refOrient) in sorted(dict_refFeatures[refId], key=lambda tup: tup[0]): 
			# Ref. annotation location
			fh_out.write (str(refId) + '\t' + str(refStart) + ':' + str(refEnd) + ':' + str(refOrient) + '\t')  

			# Ref. annotation details
			print_refAnnotDetails(dict_refFeatures[refId][(refStart, refEnd, refOrient)], fh_out)

			
			# isInterrupted directly
			if refId in dict_refAnnots_interrupted and (refStart, refEnd, refOrient) in dict_refAnnots_interrupted[refId] and len(dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)]) > 0: 
				fh_out.write("\tTrue\t")

				# Interrupted by

				ISinsertPrintStr = getISinsertStrToPrint(dict_ISannots, dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)], refStart, refEnd, refOrient, True, False, None, None, refId)
				fh_out.write(ISinsertPrintStr) 
				
				"""
				list_IStypes = getIStypes(dict_ISannots, ISstart, ISend, ISorient) 
				for IStype in list_IStypes: 
					sys.stdout.write(str(ISstart) + ':' + str(ISend) + ':' + IStype + ':' + ISorient + ';')
				"""

			else: 
				fh_out.write('\tFalse\t')
			
			
			# The left and right surrounding distance region.
			if th_surroundDist > 0: 

				##  Printing the effective surrounding dist.
				fh_out.write('\t')
				if (refStart, refEnd, refOrient) in dict_distToCheckForInsert[refId]: 
					
					[dist_left, dist_right] = dict_distToCheckForInsert[refId][(refStart, refEnd, refOrient)]
					
					fh_out.write('(') 
					if dist_left == None: 
						fh_out.write(str(dist_left))
					else: 
						fh_out.write(str(dist_left) + ":" + str(refStart) + ':seqlen=' + str(refStart - dist_left + 1))
					fh_out.write('),') 

					fh_out.write('(') 
					if dist_right == None: 
						fh_out.write(str(dist_right))
					else: 
						fh_out.write(str(refEnd) + ":" + str(dist_right) + ':seqlen=' + str(dist_right - refEnd + 1))
					fh_out.write(')') 

				## Is left interrupted
				fh_out.write('\t')
				if (refStart, refEnd, refOrient) in dict_insertsInDist[refId]: 

					if 'left' in dict_insertsInDist[refId][(refStart, refEnd, refOrient)] and len(dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['left']) > 0: 
						fh_out.write('True\t')
						ISinsertPrintStr = getISinsertStrToPrint(dict_ISannots, dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['left'], refStart, refEnd, refOrient, False, True, True, False, refId)
						fh_out.write(ISinsertPrintStr) 
						fh_out.write('\t')

					else: 
						fh_out.write('False\t\t')

					if 'right' in dict_insertsInDist[refId][(refStart, refEnd, refOrient)] and len(dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['right']) > 0:
						fh_out.write('True\t')
						ISinsertPrintStr = getISinsertStrToPrint(dict_ISannots, dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['right'], refStart, refEnd, refOrient, False, True, False, True, refId)
						fh_out.write(ISinsertPrintStr) 
						fh_out.write('\t')
					else: 
						fh_out.write('False\t\t')
			
			fh_out.write('\n')

def getIStypes(dict_ISannots, refId, ISstart, ISend, ISorient): 
	
	list_ISids = [] # [ISid1, ISid2, ...]
	list_isNew = [] # [T, F, T, ...]

	
	if refId in dict_ISannots: 
		for feature in dict_ISannots[refId][(ISstart, ISend, ISorient)]: 
			for key in feature.qualifiers: 
				# print (feature.qualifiers[key])
				if key == 'Name': 
					for aName in feature.qualifiers[key]: 
						list_ISids.append(aName)
				if key == 'isNew': 
					for aIsNew in feature.qualifiers[key]: 
						list_isNew.append(aIsNew)

					# print ('A isNew found! ')
					# print (feature.qualifiers[key])
	if len(list_ISids) != len(list_isNew):
		sys.exit('Error: Number of ISids and isNew do not match in GFF3 feature. Please check data, or code.')

	 

	return (list_ISids, list_isNew)
		


def print_refAnnotDetails(list_features, fh_out):

	isPrintColon = False

	def printColonIfTrue(): 
		if isPrintColon == True: 
			fh_out.write(':')

	for feature in list_features: 

		if feature.type: 
			fh_out.write('type=' + feature.type)
			isPrintColon = True 

		if feature.id: 

			printColonIfTrue()
			fh_out.write('id=' + feature.id) 
			isPrintColon = True

		for key in feature.qualifiers: 

			printColonIfTrue() 
			fh_out.write(key + '=') 
			for idx in range(0, len(feature.qualifiers[key])): 
				fh_out.write(feature.qualifiers[key][idx])

		fh_out.write(';')

		# print("\nSTART")
		# print (feature)
		# print("END")

def getISinsertStrToPrint(dict_ISannots, list_ISloc, refStart, refEnd, refOrient, isDirInsert, isDist, isLeft, isRight, refId): 
	printStr = '' 

	for (ISstart, ISend, ISorient) in list_ISloc: 
		(list_IStypes, list_isNew) = getIStypes(dict_ISannots, refId, ISstart, ISend, ISorient) 
		for idx_IStype, IStype in enumerate(list_IStypes): 
			
			printStr = printStr + str(ISstart) + ':' + str(ISend) + ':' + IStype + ':' + ISorient
			

			if isDirInsert == True: 

				# Is same orientation with IS.
				# if ISorient == refOrient: 
				#	printStr = printStr + ':isSameOrient=True'

				if refOrient == '+': 
					# Insert position
					insertPosition = 0
					if (ISstart >= refStart and ISstart <= refEnd):
						insertPosition = ISstart - refStart + 1 
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) ## For debugging #  + ':_1'
					elif (ISend >= refStart and ISend <= refEnd): 
						insertPosition = ISend - refStart + 1 
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) # + ':_2'
					else: 
						insertPosition = ISstart - refStart 
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) # + ':_3'
						


				elif refOrient == '-': 
					insertPosition = 0
					if (ISend <= refEnd and ISend >= refStart): 
						insertPosition = refEnd - ISend + 1
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) # + ':_4'
					elif (ISstart >= refStart and ISstart <= refEnd): 
						insertPosition = refEnd - ISstart + 1
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) # + ':_5'
					else: 
						insertPosition = refEnd - ISend
						printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ":insertPos=" + str(insertPosition) # + ':_6'
				
				 

			elif isDist == True: 
				if isLeft == True: 
					leftDist = refStart - ISend
					printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ':leftDist=' + str(leftDist) 

				if isRight == True: 
					rightDist = ISstart - refEnd
					printStr = printStr + ":isNew=" + list_isNew[idx_IStype] + ':rightDist=' + str(rightDist)
					
			printStr = printStr + ';'

	return printStr

################################################## AUX 
def listISWhichInterruptAnnots(dict_refAnnots_interrupted): 
	dict_ISdirectlyInterup = dict() # dict_[refId] => [(ISstart, ISend, ISorient)]

	for refId in dict_refAnnots_interrupted: 

		if refId not in dict_ISdirectlyInterup: 
			dict_ISdirectlyInterup[refId] = list() 

		for refInfo in dict_refAnnots_interrupted[refId]: 
			for ISinfo in dict_refAnnots_interrupted[refId][refInfo]: 
				
				if ISinfo not in dict_ISdirectlyInterup[refId]: 
					dict_ISdirectlyInterup[refId].append(ISinfo)

	return dict_ISdirectlyInterup

def getInterruptionsInDist(dict_distToCheckForInsert, dict_ISannots, dict_refAnnots_interrupted): 

	dict_IStoExclude = listISWhichInterruptAnnots(dict_refAnnots_interrupted) 

	dict_insertsInDist = dict() # dict_[refId] => [(start, end, orient)] => {left: [(ISstart, ISend, ISorient), ...], right: [(ISstart, ISend, ISorient), ...]}

	for refId in dict_distToCheckForInsert: 

		dict_insertsInDist[refId] = dict() 


		for (refStart, refEnd, refOrient) in sorted(dict_distToCheckForInsert[refId], key=lambda tup: tup[0]): 

			dict_insertsInDist[refId][(refStart, refEnd, refOrient)] = {'left': [], 'right': []}

			if refId in dict_ISannots: 
				for (start_IS, end_IS, orient_IS) in sorted(dict_ISannots[refId], key=lambda tup: tup[0]): 

					if refId in dict_IStoExclude and (start_IS, end_IS, orient_IS) not in dict_IStoExclude[refId]: 

						[dist_left, dist_right] = dict_distToCheckForInsert[refId][(refStart, refEnd, refOrient)]
						
						if dist_left != None: 
							# Check for overlap between (dist_left and refStart)

							if isOverlap(dist_left, refStart, start_IS, end_IS): 

								if (start_IS, end_IS, orient_IS) not in dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)]: 
									dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['left'].append((start_IS, end_IS, orient_IS))
								

						if dist_right != None: 
							# Check for overlap between (refEnd and dist_right) 
							
							if isOverlap(refEnd, dist_right, start_IS, end_IS): 
								if (start_IS, end_IS, orient_IS) not in dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)]: 
									dict_insertsInDist[refId][(refStart, refEnd, refOrient)]['right'].append((start_IS, end_IS, orient_IS))

			
	return dict_insertsInDist
	 


def calcEffectiveDistToCheck(dict_refFeatures, th_surroundDist): 
	
	dict_distToCheckForInsert = dict() # dict_[(start, end, orient)] = [leftDist, rightDist] #leftDist & rightDist are None if they should not be checked. 


	for refId in dict_refFeatures: 	
		dict_distToCheckForInsert[refId] = dict() 

		for (start1, end1, orient1) in sorted(dict_refFeatures[refId], key=lambda tup: tup[0]): 
			
			
			# dict_distToCheckForInsert[(start1, end1, orient1)] = [dist_left, dist_right]


			dist_left = start1 - th_surroundDist
			if dist_left < 0: 
				dist_left = 0
			dist_right = end1 + th_surroundDist

			dict_distToCheckForInsert[refId][(start1, end1, orient1)] = [(dist_left, dist_right)]

			for (start2, end2, orient2) in sorted(dict_refFeatures[refId], key=lambda tup: tup[0]): 

				
				if dist_left == None and dist_right == None: # No surrounding distance. 
					continue


				if start1 == start2 and end1 == end2 and orient1 == orient2:  # Ignore self-annotation.
					continue

				if isOverlap((start1 - th_surroundDist), (end1 + th_surroundDist), start2, end2):   
					
					if start2 >= start1 and end2 <= end1: # (s2, e2) is in the middle of (s1, e1) - effectively check everything in the surrounds
						pass 

					elif start1 >= start2 and end1 <= end2: # (s1, e1) is a sub-annotation of another annotation, no need to check the surrounds. 

						dist_right = None 
						dist_left = None 
						break


					else: 
						if dist_left != None: # Recalculate dist_left 
							if end2 >= start1 and end2 <= (end1): #  + th_surroundDist): # There is no left dist between the two. 
								dist_left = None 
							elif end2 >= dist_left and end2 <= start1:

								if (end2 + 1) > dist_left: 
									dist_left = end2 + 1
								
						if dist_right != None: # Recalculate dist_right 
							if start2 >= start1 and start2 <= end1: 
								dist_right = None
							elif start2 >= end1 and start2 <= dist_right:

								if (start2 - 1) < dist_right: 
									dist_right = start2 - 1


			if dist_left != None and dist_left >= start1: 
				dist_left = None

			if dist_right != None and dist_right <= end1: 
				dist_right = None
			

			dict_distToCheckForInsert[refId][(start1, end1, orient1)] = [dist_left, dist_right]

			#sys.stdout.write (str(start1) + ":" + str(end1) + ':' + orient1 )
			#sys.stdout.write ('\t' + 'Dist_left: ' + str(dist_left) + ' Dist_right: ' + str(dist_right) + '\n')



	#for key in dict_distToCheckForInsert: 
	#	print (str(key) + '\t' + str(dict_distToCheckForInsert[key]))

	return dict_distToCheckForInsert

def getInterruptedAnnots(dict_refFeatures, dict_ISannots): 

	dict_refAnnots_interrupted = dict() # dict[refId] => [(refStart, refEnd, refOrient)] => [(ISstart1, ISend1, ISorient1), (...), ...]
	

	for refId in dict_refFeatures: 
		for (refStart, refEnd, refOrient) in sorted(dict_refFeatures[refId], key=lambda tup: tup[0]):
			
			if refId not in dict_refAnnots_interrupted: 
				dict_refAnnots_interrupted[refId] = dict() 

			dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)] = [] 

			if refId in dict_ISannots: 
				for (ISstart, ISend, ISorient) in sorted(dict_ISannots[refId], key=lambda tup: tup[0]):    

					if isOverlap(refStart, refEnd, ISstart, ISend) == True: 
						dict_refAnnots_interrupted[refId][(refStart, refEnd, refOrient)].append((ISstart, ISend, ISorient))
						# add ISannot to interruping (so that can be excluded from further analysis)



			# print (str(refStart) + ' ' + str(refEnd) + ' ' + str(refOrient) + '\t' + str(interruptedBy)) 

	return dict_refAnnots_interrupted


def isOverlap(refStart, refEnd, ISstart, ISend):

	if (refStart >= ISstart and refStart <= ISend) or (refEnd >= ISstart and refEnd <= ISend) or (ISstart >= refStart and ISstart <= refEnd) or (ISend >= refStart and ISend <= refEnd): 
		return True 

	return False 
	 

def load_ISannot(fn_ISannotation, dict_refSeq): 

	dict_ISannots = dict() # dict_ISannots[refId] => [(loc_start, loc_end, strand)] = [gffLine, gffLine, ...]
	with open(fn_ISannotation, 'r') as fh: 
		for record in GFF.parse(fh, base_dict=dict_refSeq): 
			for feature in record.features: 

				loc_start = int(feature.location.start) + 1
				loc_end = int(feature.location.end) 

				# print (str(loc_start)  + ':' + str(loc_end))
				# print ('Strand: ' + str(feature.location.strand))
				if feature.location.strand == 1:
					strand = '+' 
				elif feature.location.strand == -1:
					strand = '-'
				else: 
					strand = '.'

				if record.id not in dict_ISannots: 
					dict_ISannots[record.id] = dict() 
				if (loc_start, loc_end, strand) not in dict_ISannots[record.id]: 
					dict_ISannots[record.id][(loc_start, loc_end, strand)] = list() 

				dict_ISannots[record.id][(loc_start, loc_end, strand)].append(feature) 


	return dict_ISannots
		


def load_refAnnot(fn_reference, fn_refAnnotation, allowedFeatTypes): 
	 
	# examiner = GFFExaminer()
	dict_refFeatures = dict() # dict_[refId] => [(start, end, orient)] => [gff3_loaded_feature1, gff3_loaded_feature2, ...]

	dict_refSeq = dict() 
	with open(fn_reference, 'r') as fh: 
		dict_refSeq = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

	with open (fn_refAnnotation, 'r') as fh: 
		for record in GFF.parse(fh , base_dict=dict_refSeq):
			for feature in record.features:
				if feature.type in allowedFeatTypes: 
					# print (str(feature.location.start) + ' ' + str(feature.location.end) + ' ' + str(feature.location.strand) + ' ' + feature.type) 

					loc_start = int(feature.location.start) + 1
					loc_end = int(feature.location.end) 

					# print(str(loc_start) + ':' + str(loc_end)) 
					strand = '+' if feature.location.strand == 1 else '-'
					if record.id not in dict_refFeatures: 
						dict_refFeatures[record.id] = dict() 
					if (loc_start, loc_end, strand) not in dict_refFeatures[record.id]: 
						dict_refFeatures [record.id][(loc_start, loc_end, strand)] = list() 
					dict_refFeatures [record.id][(loc_start, loc_end, strand)].append(feature)


					
	return (dict_refFeatures, dict_refSeq)


################################################## MAIN 
def main(): 

	
	parser = argparse.ArgumentParser(description='Identify annotations directly (hard-overlap) disrupted by IS, and within a provided threshold, IS insertion distance from the annotation. In the latter, if there is another annotation the distance is adjusted to another annotation if that comes first - so the effective distance is not the distance provided by the user, but rather the distance until the next annotation.')

	parser.add_argument('--reference', required=True, nargs=1, help='Sequence on which the annotations are defined.')
	parser.add_argument('--ref_annotation', required=True, nargs=1, help='Annotations in the .gff3 format. This file can be downloaded from NCBI. You should not provide the line containing the annotation of the reference itself (e.g by removing it from the gff3 file).')
	parser.add_argument('--IS_annotation', required=True, nargs=1, help='IS annotations generated by WaIS in the .gff3 format.')

	parser.add_argument('--surrounding_distance', required=False, nargs=1, default=[100], help='The distance (in base pairs) around which the IS insertions should be reported (such as for insertions in promotor regions). To switch this feature off, set this value to 0. Default=100.')
	parser.add_argument('--annotation_types', help='Types of features to calculate IS insertions in. Default=[\'gene\', \'pseudogene\', \'riboswitch\'].', default=['gene', 'pseudogene', 'riboswitch'], nargs='+')

	parser.add_argument('--outfile', required=True, nargs=1, help="Name of output file to save disrupted insertions to.")
	parser.add_argument('--outfile_onlyInterrupted', required=False, default=[None], nargs=1, help='Output in alternative format (printing only those interrupted, & switched left and right to upstream and downstream.')
	
	args = parser.parse_args()

	sys.stderr.write("Note: For a gene where isInterrupted==True, insertion at the first position of the gene is recorded as 1 (i.e. it is 1 indexed), and one position before the gene is recorded as -1. There is never a 0.\n\n")
	# print (args.annot_types)

	insertionsWrtRefAnnotations(args.reference[0], args.ref_annotation[0], args.IS_annotation[0], args.annotation_types, int(args.surrounding_distance[0]), args.outfile[0], args.outfile_onlyInterrupted[0]) 

if __name__=='__main__': 
	main() 