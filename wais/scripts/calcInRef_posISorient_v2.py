#!/usr/bin/python3

# from Scripts.mergeLocalCounts import SIDE_L, SIDE_S
import calcInContig_posISOrient as calcInContig
import re
import argparse
import sys 
import calcInContig_posISOrient
from collections import OrderedDict

Col_qSeqId = 0; 
Col_sSeqId = 1; 
Col_alignLen = 3
Col_pIdent = 2; 
Col_qSeqLen = 12 
Col_sSeqLen = 13

Col_qStart = 6
Col_qEnd = 7
Col_sStart = 8 
Col_sEnd = 9

Col_gff3_info = 8

SITES = 'sites'
INFO = 'information'
ISinContigsIds = 'ISids_inContigs'

DIRECT_TO_REF = 'DirectInRef'
# CALC_FROM_CONTIG = 'CalcFromContig'
DIRECT_TO_CONTIG = 'DirectInContig'
READ_MAPPING = 'ReadMapping'

EDGE = 'Edge'
MIDDLE = 'Middle'
SIDE_L = 'Side_l'
SIDE_S = 'Side_s'

#################################### TOP_LVL 
def calcPosToRef(fn_contigToRef, fn_IStoContig, th_alignDiff_IStoContig, fn_gff_ISannot, fn_directIStoRef, fn_out, isMerged, th_distForMerge, isIgnoreOrient, isIgnoreIStype, separator, th_directOverlap, th_maxAlignLenDiff, th_contigToRef_ISonlyLenDiff):

	# 1. Load direct IStoRef
	dict_direct = load_direct_IStoRef(fn_directIStoRef, isMerged, isIgnoreOrient, isIgnoreIStype, separator)
	
	
	# """
	for (refId, refLen) in dict_direct: 
		for (ISid, refStart, refEnd, orient, isContigFlipped) in dict_direct[(refId, refLen)]: 
			print (refId + '\t' + str(refLen) + '\t' + ISid + '\t' + str(refStart) + '\t' + str(refEnd) + '\t' + orient + '\t'  + str(isContigFlipped))
	# """

	# return 

	# 2. contigs which are completely IS 
	list_contigsWithComplIS = identifyISonlyContigs(fn_IStoContig, th_alignDiff_IStoContig)


	# 4. Load IS in Contig (gff)
	(dict_IStoContig, dict_color) = load_ISinContig_gff(fn_gff_ISannot)


	# 3. load contigToRef 
	(dict_contigToRef, list_allContigs, dict_contigToRef_alignLen) = load_contigToRef(fn_contigToRef, list_contigsWithComplIS, dict_direct, th_maxAlignLenDiff, dict_IStoContig, th_contigToRef_ISonlyLenDiff) 

	 
	# 5. Print calc and print IS position from contig to ref
	list_ISinRef_calc = determineISforRef(dict_contigToRef, dict_IStoContig)

	
	if isMerged == 'True':
		(dict_mergedIgnoreIStype, dict_mergedIgnoreOrient, dict_mergedIStoRef, list_seqRegions) = doTheMerge(dict_direct, list_ISinRef_calc, th_distForMerge, isIgnoreOrient, isIgnoreIStype, separator)
		
		if len(dict_mergedIgnoreIStype) > 0: 
			print('Printing for --ignoreIStype')
			printGff3Merged(dict_mergedIgnoreIStype, fn_out, dict_direct, dict_color, list_seqRegions, True, True, th_directOverlap)
			# getISnotMappedToRef()
		elif len(dict_mergedIgnoreOrient) > 0:
			print('Printing for --ignoreOrient')
			printGff3Merged(dict_mergedIgnoreOrient, fn_out, dict_direct, dict_color, list_seqRegions, True, False, th_directOverlap)
		else: 
			print('Printing for --merged')
			dict_ISidsInContigs_counts = printGff3Merged(dict_mergedIStoRef, fn_out, dict_direct, dict_color, list_seqRegions, False, False, th_directOverlap)

			print ('dict_ISidsInContigs_counts is .... ')
			for uniqId in dict_ISidsInContigs_counts: 
				print (uniqId + '\t' + str(dict_ISidsInContigs_counts[uniqId]))
			# print (dict_ISidsInContigs_counts)
			# getISnotMappedToRef(list_allContigs, dict_contigToRef_alignLen, dict_IStoContig, list_contigsWithComplIS)

			# calcAndPrintISnotMappedToRef(dict_ISidsInContigs_counts, dict_IStoContig)
		# Printing the 

	# if isMerged == 'False': 
	else: 

		fh_out = open(fn_out, 'w+')
		list_refsEncountered = [] 
		
		for (refId, refLen) in dict_direct: 

			if ((refId, refLen) not in list_refsEncountered): 
				list_refsEncountered.append((refId, str(refLen)))


		

			for (ISid, start_ISinRef_dir, end_ISinRef_dir, orient_IStoRef_dir, isRefFlipped) in dict_direct[(refId, refLen)]: 

				if ISid not in dict_color:
					dict_color[ISid] = calcInContig_posISOrient.generateRandomColor()

					 
				# calcInContig.printGff3Line(fh_out, refId, calcInContig.PROG_NAME, calcInContig.SO_, start_ISinRef_dir, end_ISinRef_dir, 'Inf', orient_IStoRef_dir, '.', 'Name=' + ISid + ';isCalcFromContig=False' + ';color=' + dict_color[ISid]); 



		for (refInfo, progName, seqFeatType, calcStart_IStoRef, calcEnd_IStoRef, score, calcOrient_IStoRef, phase, attributes)  in list_ISinRef_calc: 
			# if not isInDirect(dict_direct, int(calcStart_IStoRef), int(calcEnd_IStoRef)):  

			if (refInfo not in list_refsEncountered): 
				list_refsEncountered.append(refInfo)


			## TODO: Check here ifInDirect (?)
 			# calcIfInDirect(dict_direct, refInfo, calcStart_IStoRef, calcEnd_IStoRef, calcOrient_IStoRef, attributes)

			calcInContig.printGff3Line(fh_out, refInfo[0], progName, seqFeatType, calcStart_IStoRef, calcEnd_IStoRef, score, calcOrient_IStoRef, phase, attributes)

		
		fh_out.close() 



		with open(fn_out, 'r') as original: 
			data = original.read()
		
		with open(fn_out, 'w') as modified: 
			modified.write('##gff-version 3.1.26' + '\n')
			for (refId, refLen) in list_refsEncountered: 
				modified.write('##sequence-region ' + refId + ' ' + '1' + ' ' + str(refLen) + '\n')
			
			modified.write(data)

		modified.close() 
	
	# for key in list_ISinRef_calc: 
	#    print (str(key))
	"""
	for key in dict_contigToRef: 
		for anAlign in dict_contigToRef[key]: 
			print (str(key) + "\t" + str(anAlign))
	"""
	
	"""
	for key in dict_direct: 
		for anAlign in dict_direct[key]: 
			print (str(key) + "\t" + str(anAlign)) 
	"""


#################################### AUX
def checkIfMerge_directRefIS(list_refInsertions, isMerged, isIgnoreOrient, isIgnoreIStype, separator, ISid_1, refStart_1, refEnd_1, orient_1, isContigFlipped_1):

	isUpdated = False
	for idx, (ISid_2, refStart_2, refEnd_2, orient_2, isContigFlipped_2) in enumerate(list_refInsertions):

		if isIgnoreIStype == 'True' and isAnyOverlap(refStart_1, refEnd_1, refStart_2, refEnd_2, 0): 
			# get the positions, merge the ISids and orient
			(s_new, e_new) = getNewPosForMerge_directRefIS(refStart_1, refEnd_1, refStart_2, refEnd_2)
			orient_new = orient_1 if orient_1 == orient_2 else '.'
			ISid_new = getNewISid_directRefIS(ISid_1, ISid_2, separator)
			list_refInsertions[idx] = (ISid_new, s_new, e_new, orient_new, isContigFlipped_1) 

			isUpdated = True 
			break 

		elif isIgnoreOrient == 'True' and ISid_1 == ISid_2 and isAnyOverlap(refStart_1, refEnd_1, refStart_2, refEnd_2, 0): 

			# Merge the orient and the positions 
			(s_new, e_new) = getNewPosForMerge_directRefIS(refStart_1, refEnd_1, refStart_2, refEnd_2)
			
			orient_new = orient_1 if orient_1 == orient_2 else '.'
			list_refInsertions[idx] = (ISid_2, s_new, e_new, orient_new, isContigFlipped_1) 

			isUpdated = True 
			break 

		elif isMerged == 'True' and ISid_1 == ISid_2 and orient_1 == orient_2 and isAnyOverlap(refStart_1, refEnd_1, refStart_2, refEnd_2, 0): 
		# Merge the positions

			# print ('Do the merge!!!')
			(s_new, e_new) = getNewPosForMerge_directRefIS(refStart_1, refEnd_1, refStart_2, refEnd_2)
			list_refInsertions[idx] = (ISid_2, s_new, e_new, orient_1, isContigFlipped_1) 

			isUpdated = True 
			break 

	return isUpdated

def getNewPosForMerge_directRefIS(s1, e1, s2, e2):

	s_new = s1 if s1 <= s2 else s2
	e_new = e1 if e1 >= e2 else e2 

	return (s_new, e_new)   


def isISonlyAlignedRegion(contigId, contigStart, contigEnd, dict_IStoContig, th_contigToRef_ISonlyLenDiff): 

	Gff3_start = 3
	Gff3_end = 4
	Gff3_orient = 5

	alignLen = contigEnd - contigStart 


	# print ('############### IS_TO_CONTIG')
	if contigId in dict_IStoContig: 
		for arr_anIns in dict_IStoContig[contigId]: 
			
			anIns_start = int(arr_anIns[Gff3_start]) 
			anIns_end = int(arr_anIns[Gff3_end])

			isOverlap = isAnyOverlap(anIns_start, anIns_end, contigStart, contigEnd, 0) # hard-overlap.
			
			anIns_alignLen = anIns_end - anIns_start 
			
			if isOverlap == True and abs(alignLen - anIns_alignLen) <= th_contigToRef_ISonlyLenDiff: 
				# print (str(anIns_start) + ':' + str(anIns_end) + '\t' + str(contigStart) + ':' + str(contigEnd) + '\t' + str(anIns_alignLen) + '\t' + str(alignLen) + '\t' + str(abs(alignLen - anIns_alignLen)))

				return True

			# print ('ContigId found in dict_ISinContig ' + contigId)

	"""	
	for contigId in dict_IStoContig: 
		print (contigId) 
		for arr_anIns in dict_IStoContig[contigId]: 
			
			print (arr_anIns[Gff3])
	"""
	return False 


def extractISidFromAttr(attrs): 

	arr = attrs.split(";")
	ISid = '';
	directContig = 0; 
	readMapping = 0;
	edge = 0
	middle = 0; 
	side_l = 0
	side_s = 0
	uniqId = -1; 

	for val in arr: 
		if re.match(r'^Name\=', val): 
			ISid = re.sub('Name=', '', val)

		elif re.match(r'^num_Inf\=', val): 
			directContig = re.sub('num_Inf=', '', val)
			directContig = int(directContig)

		elif re.match(r'^combined_Depth\=', val): 
			readMapping = re.sub('combined_Depth=', '', val)
			readMapping = int(readMapping)

		elif re.match(r'^side_s\=', val): 
			side_s = re.sub('side_s=', '', val)
			side_s = int(side_s)

		elif re.match(r'^side_l\=', val): 
			side_l = re.sub('side_l=', '', val)
			side_l = int(side_l)

		elif re.match(r'^positionInContig\=', val): 
			posInContig = re.sub('positionInContig=', '', val)
			if (posInContig == EDGE): 
				edge = 1 
			elif posInContig == MIDDLE: 
				middle = 1 

		elif re.match(r'^uniqId\=', val):
			uniqId = re.sub(r'^uniqId\=', '', val)
			# print ('Found uniqId as ' + val)




		""" 
		if re.match('^Name\=', val): 
			ISid = re.sub('Name=', '', val)

		if re.match('^Name\=', val): 
			ISid = re.sub('Name=', '', val)
		"""

	return (ISid, directContig, readMapping, edge, middle, side_s, side_l, uniqId)

def findIfInDirect(dict_direct, start_calc, end_calc, isIgnoreOrient, isIgnoreIStype, ISid_calc, orient_calc, th_directOverlap):

	in_refStart = None
	in_refEnd = None
	in_ISid = None
	in_orient = None

	isInDirect = False 

	list_refIS = []

	for (refId, refLen) in dict_direct: 
		for (ISid, start_ISinRef_dir, end_ISinRef_dir, orient_IStoRef_dir, isRefFlipped) in dict_direct[(refId, refLen)]:  
			if ISid_calc == ISid and orient_calc == orient_IStoRef_dir and isAnyOverlap(start_calc, end_calc, start_ISinRef_dir, end_ISinRef_dir, th_directOverlap): 

				isInDirect = True 

				in_refStart = start_ISinRef_dir
				in_refEnd = end_ISinRef_dir
				in_ISid = ISid
				in_orient = orient_IStoRef_dir

				# return (True, in_refStart, in_refEnd, in_ISid, in_orient)
				list_refIS.append((in_refStart, in_refEnd, in_ISid, in_orient))


			elif isIgnoreIStype == True and isAnyOverlap(start_calc, end_calc, start_ISinRef_dir, end_ISinRef_dir, th_directOverlap): 
				
				isInDirect = True 

				in_refStart = start_ISinRef_dir
				in_refEnd = end_ISinRef_dir
				in_ISid = ISid
				in_orient = orient_IStoRef_dir

				# return (True, in_refStart, in_refEnd, in_ISid, in_orient)
				list_refIS.append((in_refStart, in_refEnd, in_ISid, in_orient))
			
			
			elif isIgnoreOrient == True and ISid_calc == ISid and isAnyOverlap(start_calc, end_calc, start_ISinRef_dir, end_ISinRef_dir, th_directOverlap): 
				
				isInDirect = True 

				in_refStart = start_ISinRef_dir
				in_refEnd = end_ISinRef_dir
				in_ISid = ISid
				in_orient = orient_IStoRef_dir

				# return (True, in_refStart, in_refEnd, in_ISid, in_orient) 
				list_refIS.append((in_refStart, in_refEnd, in_ISid, in_orient))


			


	return (isInDirect, list_refIS)



def doTheMerge(dict_direct, list_ISinRef_calc, th_distForMerge, isIgnoreOrient, isIgnoreIStype, separator):

	list_seqRegions = [] 

	dict_mergedIStoRef = dict() # dict[(refId)] => [(start, end, ISid, orient)] => [(...)]
	dict_mergedIgnoreOrient = dict() 
	dict_mergedIgnoreIStype = dict()

	
	for (refId, refLen) in dict_direct:
		if (refId, int(refLen)) not in list_seqRegions: 
			list_seqRegions.append((refId, int(refLen))) 
			
	"""
		for (ISid, start_ISinRef_dir, end_ISinRef_dir, orient_IStoRef_dir, isRefFlipped) in sorted(dict_direct[(refId, refLen)],  key=lambda tup: tup[1]): 

			checkAndAddToMerge(dict_mergedIStoRef, refId, start_ISinRef_dir, end_ISinRef_dir, ISid, orient_IStoRef_dir, th_distForMerge, 1, 0, 0, 0, 0, 0, 0)
	"""

	for ((refId, refLen), gff3_tool, gff3_IS, refStart, refEnd, score, orient, gff_col, attrs) in sorted(list_ISinRef_calc, key=lambda tup: tup[3]): 
		
		if (refId, int(refLen)) not in list_seqRegions: 
			list_seqRegions.append((refId, int(refLen))) 

			
		(ISid, directContig, readMapping, edge, middle, side_s, side_l, uniqId) = extractISidFromAttr(attrs)
		
		# print(attrs)
		# print (str((ISid, directContig, readMapping, edge, middle, side_s, side_l)))
		# print ('Val: ' + str(refStart) + ' ' + str(refEnd) + ' ' + str(orient) + ' ' + ISid)
		# print(refInfo)
		checkAndAddToMerge(dict_mergedIStoRef, refId, refStart, refEnd, ISid, orient, th_distForMerge, 0, directContig, readMapping, edge, middle, side_s, side_l, uniqId)

		

	if isIgnoreOrient == 'True' or isIgnoreIStype == 'True': 
		dict_mergedIgnoreOrient = reduceSet_3(dict_mergedIStoRef, th_distForMerge, separator) 

	if isIgnoreIStype == 'True': 
		dict_mergedIgnoreIStype = reduceSet_4(dict_mergedIgnoreOrient, th_distForMerge, separator)


	print ("################################")

	
	if len(dict_mergedIgnoreIStype) > 0: 
		siteKeys = []
		print("Printing dict_mergedIgnoreIStype: ")
		for key in dict_mergedIgnoreIStype: 
			for site in sorted(dict_mergedIgnoreIStype[key][SITES]): 
				print ("Key: " + str(site))
				# print (dict_mergedIStoRef[key][INFO].keys())
				# if (site[0], site[1]) in dict_mergedIgnoreOrient[key][INFO]: 
				print(dict_mergedIgnoreIStype[key][INFO][site])
				# siteKeys.append((site[0], site[1]))

		print ("Total keys: " + str(len(dict_mergedIgnoreIStype[key][SITES])) + ' ' + str(len(dict_mergedIgnoreIStype[key][INFO])))

		# diffKeys = list(set(dict_mergedIgnoreIStype[key][INFO].keys()) - set(siteKeys))

		# print ('The diffKeys are: ')
		# print (diffKeys)
			
	elif len(dict_mergedIgnoreOrient) > 0:
		siteKeys = []
		print("Printing dict_mergedIgnoreOrient: ")
		for key in dict_mergedIgnoreOrient: 
			for site in sorted(dict_mergedIgnoreOrient[key][SITES]): 
				print ("Key: " + str(site))
				# print (dict_mergedIStoRef[key][INFO].keys())
				# if (site[0], site[1]) in dict_mergedIgnoreOrient[key][INFO]: 
				print(dict_mergedIgnoreOrient[key][INFO][site])
				# siteKeys.append((site[0], site[1]))

		print ("Total keys: " + str(len(dict_mergedIgnoreOrient[key][SITES])) + ' ' +  str(len(dict_mergedIgnoreOrient[key][INFO].keys())))

		diffKeys = list(set(dict_mergedIgnoreOrient[key][INFO].keys()) - set(dict_mergedIgnoreOrient[key][SITES]))

		print ('The diffKeys are: ')
		print (diffKeys)

	else: 
	
		siteKeys = []
		for key in dict_mergedIStoRef: 
			for site in sorted(dict_mergedIStoRef[key][SITES]): 
				print ("Key: " + str(site))
				# print (dict_mergedIStoRef[key][INFO].keys())
				# if site in dict_mergedIStoRef[key][INFO]: 
				print(dict_mergedIStoRef[key][INFO][site])

				# siteKeys.append()

		print ("Total keys: " + str(len(dict_mergedIStoRef[key][SITES])) + ' ' + str(len(dict_mergedIStoRef[key][INFO])))

		# diffKeys_1 = list(set(dict_mergedIStoRef[key][INFO].keys()) - set(dict_mergedIStoRef[key][SITES]))

		# diffKeys_2 = list(set(dict_mergedIStoRef[key][SITES]) - set(dict_mergedIStoRef[key][INFO].keys()))


		# print ('The diffKeys are: ')
		# print (diffKeys_1)
		# print (diffKeys_2)


	return (dict_mergedIgnoreIStype, dict_mergedIgnoreOrient, dict_mergedIStoRef, list_seqRegions)

def extractISidAndColor(arrLine, dict_colors):

	# print ("The arrLine is ") 
	# print (arrLine[Col_gff3_info]) 

	arr = arrLine[Col_gff3_info].split(';')
	ISid = '' 
	color = '' 
	for val in arr: 
		if 'Name=' in val: 
			ISid = re.sub('^Name=', '', val)

		if 'color=' in val: 
			color = re.sub('^color=', '', val)


	if ISid not in dict_colors:
		if color != '':  
			dict_colors[ISid] = color
		else: 
			## 
			dict_colors[ISid] = calcInContig_posISOrient.generateRandomColor()

	return (ISid, color)

def printGff3Merged(dict_, fn_out, dict_direct, dict_color, list_seqRegions, isIgnoreOrient, isIgnoreIStype, th_directOverlap): 

	dict_ISidsInContigs_counts = dict() 


	fh_out = open (fn_out, 'w+')
	
	fh_out.write('##gff-version 3.1.26' + '\n')
	
	for (refId, refLen) in list_seqRegions: 
		if refId in dict_: 
			fh_out.write('##sequence-region ' + refId + ' ' + '1' + ' ' + str(refLen) + '\n')

	for refId in dict_: 
		for (startPos, endPos, ISid, orient) in dict_[refId][SITES]: 

			dict_info = dict_[refId][INFO][(startPos, endPos, ISid, orient)]
			
			# print ('The dict_info is')
			# print (dict_info)

			attributes = 'Name=' + ISid + ';'
			attributes = attributes + DIRECT_TO_REF + '=' + str(dict_info[DIRECT_TO_REF]) + ';'
			attributes = attributes + DIRECT_TO_CONTIG + '=' + str(dict_info[DIRECT_TO_CONTIG]) + ';'
			attributes = attributes + READ_MAPPING + '=' + str(dict_info[READ_MAPPING]) + ';'
			attributes = attributes + EDGE + '=' + str(dict_info[EDGE]) + ';'
			attributes = attributes + MIDDLE + '=' + str(dict_info[MIDDLE]) + ';'
			attributes = attributes + SIDE_L + '=' + str(dict_info[SIDE_L]) + ';'
			attributes = attributes + SIDE_S + '=' + str(dict_info[SIDE_S]) + ';'
			attributes = attributes + ISinContigsIds + '=' + ','.join(dict_info[ISinContigsIds]) + ';'

			addCounts_ISidsInContigs(dict_ISidsInContigs_counts, dict_info[ISinContigsIds])

			if (ISid not in dict_color):
				color = calcInContig_posISOrient.generateRandomColor()
				dict_color[ISid] = color

			attributes = attributes + 'color' + '=' + dict_color[ISid] + ';'

			(isInDirect, list_directRefIS) = findIfInDirect(dict_direct, startPos, endPos, isIgnoreOrient, isIgnoreIStype, ISid, orient, th_directOverlap)
			attributes = attributes + 'isNew=' + str(not isInDirect) + ';'

			if isInDirect == True: 
				for idx, (direct_start, direct_end, direct_ISid, direct_orient) in enumerate(list_directRefIS): 
					if idx == 0: 
						attributes = attributes + 'refIS='
					
						
					attributes = attributes + str(direct_start) + ':' + str(direct_end) + ':' + direct_ISid + ':' + direct_orient

					if idx != len(list_directRefIS) - 1: 
						attributes = attributes + ','

				attributes = attributes + ';'
		

			if '+' in orient and '-' in orient: 
				calcInContig.printGff3Line(fh_out, refId, 'WiIS', 'insertion_sequence', startPos, endPos, '.', '.', '.', attributes)
			else: 
				calcInContig.printGff3Line(fh_out, refId, 'WiIS', 'insertion_sequence', startPos, endPos, '.', orient, '.', attributes)

			


	fh_out.close() 

	return dict_ISidsInContigs_counts

def reduceSet_4(dict_reducedSet1, th_forMergingOverlaps, separator):

	dict_reducedSet2 = dict() # dict_reducedSet2[refId|contigId] = [(startPos|min, endPos|max, ISid, 'orient1:orient2')]

	for key in dict_reducedSet1:
		for (start_set1, end_set1, ISid_set1, orient_set1) in dict_reducedSet1[key][SITES]:

			if key not in dict_reducedSet2:
				dict_reducedSet2[key] =  {SITES: [], INFO: {}}
				dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
				#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None)

				# updateMergedSiteInfo(dict_reducedSet2[key][INFO], None, start_set1, end_set1, None, None, None, None, dict_reducedSet1[key][INFO][(start_set1, end_set1)][DIRECT_TO_REF], dict_reducedSet1[key][INFO][(start_set1, end_set1)][DIRECT_TO_CONTIG], dict_reducedSet1[key][INFO][(start_set1, end_set1)][READ_MAPPING], dict_reducedSet1[key][INFO][(start_set1, end_set1)][EDGE], dict_reducedSet1[key][INFO][(start_set1, end_set1)][MIDDLE], dict_reducedSet1[key][INFO][(start_set1, end_set1)][SIDE_L], dict_reducedSet1[key][INFO][(start_set1, end_set1)][SIDE_S]) 

				updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			else: # Loop-through, check - merge or create new

				isOverlapFound = False
				for idx, (start_set2, end_set2, ISid_set2, orient_set2) in enumerate(dict_reducedSet2[key][SITES]):
					if isAnyOverlap(start_set1, end_set1, start_set2, end_set2, th_forMergingOverlaps):

						# dict_redu
						# print ('Do the merge!')

						newStart = start_set2 if (start_set2 <= start_set1) else start_set1
						newEnd = end_set2 if end_set2 >= end_set1 else end_set1


						
						orient = separator.join(list((set(orient_set1.split(separator) + orient_set2.split(separator)))))

						ISid = separator.join(sorted(list((set(ISid_set1.split(separator) + ISid_set2.split(separator))))))	


						dict_reducedSet2[key][SITES][idx] = (newStart, newEnd, ISid, orient)
						#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2)

						updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2, ISid, orient, ISid_set1, orient_set1, ISid_set2, orient_set2)

						isOverlapFound = True
						break

				if isOverlapFound == False:
					dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
					#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None)

					updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			# 	if isAnyOverlapß
			# print ('\t' + str())

# 	pass

	return dict_reducedSet2

def reduceSet_3(dict_reducedSet1, th_forMergingOverlaps, separator):

	dict_reducedSet2 = dict() # dict_reducedSet2[refId|contigId] = {'sites' => [(startPos|min, endPos|max, ISid, 'orient1:orient2')]; 'information' => {numInf: '0', combinedScore: '', side_s: '', side_l: '', ...}]

	for key in dict_reducedSet1:
		for (start_set1, end_set1, ISid_set1, orient_set1) in sorted(dict_reducedSet1[key][SITES], key=lambda tup:tup[0]):

			if key not in dict_reducedSet2:
				dict_reducedSet2[key] = {SITES: [], INFO: {}}
				dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))

				print ("The key is " + str((start_set1, end_set1, ISid_set1, orient_set1)))
				#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None)
				# keyLocal = (start_set1, end_set1, ISid_set1, orient_set1)
				
				# updateMergedSiteInfo(dict_reducedSet2[key][INFO], None, start_set1, end_set1, None, None, None, None, dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_REF], dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_CONTIG], dict_reducedSet1[key][INFO][keyLocal][READ_MAPPING], dict_reducedSet1[key][INFO][keyLocal][EDGE], dict_reducedSet1[key][INFO][keyLocal][MIDDLE], dict_reducedSet1[key][INFO][keyLocal][SIDE_L], dict_reducedSet1[key][INFO][keyLocal][SIDE_S], ISid_set1, orient_set1, ISid_set1, orient_set1) 


				updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			else: # Loop-through, check - merge or create new

				isOverlapFound = False
				for idx, (start_set2, end_set2, ISid_set2, orient_set2) in enumerate(dict_reducedSet2[key][SITES]):
					if ISid_set1 == ISid_set2 and isAnyOverlap(start_set1, end_set1, start_set2, end_set2, th_forMergingOverlaps):

						# dict_redu
						# print ('Do the merge!')

						newStart = start_set2 if (start_set2 <= start_set1) else start_set1
						newEnd = end_set2 if end_set2 >= end_set1 else end_set1


						orient = separator.join(list((set(orient_set1.split(separator) + orient_set2.split(separator)))))

						# if (orient_set1 not in orient_set2.split(separator)):
						# 	orient = orient_set2 + separator + orient_set1

						dict_reducedSet2[key][SITES][idx] = (newStart, newEnd, ISid_set2, orient)
						#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2)

						# keyLocal = (start_set2, end_set2, ISid_set2, orient_set2)
						# updateMergedSiteInfo(dict_reducedSet2[key][INFO], None, newStart, newEnd, start_set2, end_set2, None, None, dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_REF], dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_CONTIG], dict_reducedSet1[key][INFO][keyLocal][READ_MAPPING], dict_reducedSet1[key][INFO][keyLocal][EDGE], dict_reducedSet1[key][INFO][keyLocal][MIDDLE], dict_reducedSet1[key][INFO][keyLocal][SIDE_L], dict_reducedSet1[key][INFO][keyLocal][SIDE_S], ISid_set1, orient, ISid_set2, orient_set2) 

						updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2, ISid_set1, orient, ISid_set1, orient_set1, ISid_set2, orient_set2)

						isOverlapFound = True
						break

				if isOverlapFound == False:
					dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
					#addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None)

					# keyLocal = (start_set1, end_set1, ISid_set1, orient_set1)
					# updateMergedSiteInfo(dict_reducedSet2[key][INFO], None, start_set1, end_set1, None, None, None, None, dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_REF], dict_reducedSet1[key][INFO][keyLocal][DIRECT_TO_CONTIG], dict_reducedSet1[key][INFO][keyLocal][READ_MAPPING], dict_reducedSet1[key][INFO][keyLocal][EDGE], dict_reducedSet1[key][INFO][keyLocal][MIDDLE], dict_reducedSet1[key][INFO][keyLocal][SIDE_L], dict_reducedSet1[key][INFO][keyLocal][SIDE_S], ISid_set1, orient_set1, ISid_set1, orient_set1) 


					updateMergedSiteInfo_2(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			# 	if isAnyOverlapß
			# print ('\t' + str())

# 	pass

	return dict_reducedSet2

def checkAndAddToMerge(dict_refWithMerges, refId, startOrig, endOrig, ISidOrig, orientOrig, th_distForMerge, isCalcFromContig, directContig, readMapping, edge, middle, side_s, side_l, uniqId):

	if (startOrig > endOrig): 
		tmp = startOrig
		startOrig = endOrig 
		endOrig = tmp 

		tmp = side_s 
		side_s = side_l 
		side_l = tmp 

		if orientOrig == '+': 
			orientOrig = '-'
		else: 
			orientOrig = '+'

	if refId not in dict_refWithMerges: 
		dict_refWithMerges[refId] = {SITES: [], INFO: {}}

	
		dict_refWithMerges[refId][SITES].append((startOrig, endOrig, ISidOrig, orientOrig))
		
		
		
		## updateMergedSiteInfo(dict_refWithMerges[refId][INFO], None, startOrig, endOrig, None, None, None, None, isCalcFromContig, directContig, readMapping, edge, middle, side_s, side_l, ISidOrig, orientOrig, ISidOrig, orientOrig)
		
		dict_prev = dict() 
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)] = dict()
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_REF] = isCalcFromContig
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_CONTIG] = directContig
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][READ_MAPPING] = readMapping
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][EDGE] = edge
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][MIDDLE] = middle
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_S] = side_s
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_L] = side_l
		dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][ISinContigsIds] = [uniqId]
		
		updateMergedSiteInfo_2(dict_refWithMerges[refId][INFO], dict_prev, startOrig, endOrig, startOrig, endOrig, None, None, ISidOrig, orientOrig, ISidOrig, orientOrig, None, None)

	else: 

		isOverlapFound = False
		for idx, (start_saved, end_saved, IS_saved, orient_saved) in enumerate(dict_refWithMerges[refId][SITES]): 

			# print ('startOrig:' + str(startOrig) + ' endOrig:' + str(endOrig) + ' start_saved: ' + str(start_saved) + ' end_saved: ' + str(end_saved))
			if ISidOrig == IS_saved and orientOrig == orient_saved and isAnyOverlap(startOrig, endOrig, start_saved, end_saved, int(th_distForMerge)):
				

				newStart = startOrig if (startOrig <= start_saved) else start_saved
				newEnd = endOrig if endOrig >= end_saved else end_saved 

				dict_refWithMerges[refId][SITES][idx] = (newStart, newEnd, IS_saved, orient_saved)

				
				
				
				### updateMergedSiteInfo(dict_refWithMerges[refId][INFO], None, newStart, newEnd, start_saved, end_saved, None, None, isCalcFromContig, directContig, readMapping, edge, middle, side_s, side_l, ISidOrig, orientOrig, IS_saved, orient_saved)

				dict_prev = dict() 
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)] = dict()
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_REF] = isCalcFromContig
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_CONTIG] = directContig
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][READ_MAPPING] = readMapping
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][EDGE] = edge
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][MIDDLE] = middle
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_S] = side_s
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_L] = side_l
				dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][ISinContigsIds] = [uniqId]
				updateMergedSiteInfo_2(dict_refWithMerges[refId][INFO], dict_prev, newStart, newEnd, startOrig, endOrig, start_saved, end_saved, ISidOrig, orientOrig, ISidOrig, orientOrig, ISidOrig, orientOrig)

				# print ("Updated: " + str(newStart) + ":" + str(newEnd) + ":" + str(ISidOrig) + ':" + str(or'  ' PrevKey: ' + str(start_saved) + ":" + str(end_saved) )

				isOverlapFound = True 
				break 

		if isOverlapFound == False: 
			# Add as new Site
			dict_refWithMerges[refId][SITES].append((startOrig, endOrig, ISidOrig, orientOrig))
			## updateMergedSiteInfo(dict_refWithMerges[refId][INFO], None, startOrig, endOrig, None, None, None, None, isCalcFromContig, directContig, readMapping, edge, middle, side_s, side_l, ISidOrig, orientOrig, ISidOrig, orientOrig)

			dict_prev = dict() 
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)] = dict()
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_REF] = isCalcFromContig
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][DIRECT_TO_CONTIG] = directContig
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][READ_MAPPING] = readMapping
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][EDGE] = edge
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][MIDDLE] = middle
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_S] = side_s
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][SIDE_L] = side_l
			dict_prev[(startOrig, endOrig, ISidOrig, orientOrig)][ISinContigsIds] = [uniqId]
			updateMergedSiteInfo_2(dict_refWithMerges[refId][INFO], dict_prev, startOrig, endOrig, startOrig, endOrig, None, None, ISidOrig, orientOrig, ISidOrig, orientOrig, None, None)

def updateMergedSiteInfo_2(dict_info_new, dict_info_prev, newStart, newEnd, prevStart, prevEnd, oldStart, oldEnd, newISid, newOrient, prevISid, prevOrient, oldISid, oldOrient):

	if (newStart, newEnd, newISid, newOrient) not in dict_info_new:
		dict_info_new[(newStart, newEnd, newISid, newOrient)] = {}
		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_REF] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][DIRECT_TO_REF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_CONTIG] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][DIRECT_TO_CONTIG]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][READ_MAPPING] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][READ_MAPPING]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][EDGE] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][EDGE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][MIDDLE] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][MIDDLE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_L]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][ISinContigsIds] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][ISinContigsIds]
        
	else: 
		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_REF] = dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_REF] +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][DIRECT_TO_REF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_CONTIG] = dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_CONTIG]  +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][DIRECT_TO_CONTIG]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][READ_MAPPING] = dict_info_new[(newStart, newEnd, newISid, newOrient)][READ_MAPPING]  +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][READ_MAPPING]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][EDGE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][EDGE]  +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][EDGE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][MIDDLE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][MIDDLE]  +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][MIDDLE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] + dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_L]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][ISinContigsIds] = dict_info_new[(newStart, newEnd, newISid, newOrient)][ISinContigsIds] + dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][ISinContigsIds]

	if oldStart != None and oldEnd != None and oldISid != None and oldOrient != None and not (newStart == oldStart and newEnd == oldEnd and newISid == oldISid and oldOrient == newOrient): 
		# print ("Now also add, this info, then del key.")

		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_REF] = dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_REF] +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][DIRECT_TO_REF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_CONTIG] = dict_info_new[(newStart, newEnd, newISid, newOrient)][DIRECT_TO_CONTIG]  +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][DIRECT_TO_CONTIG]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][READ_MAPPING] = dict_info_new[(newStart, newEnd, newISid, newOrient)][READ_MAPPING]  +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][READ_MAPPING]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][EDGE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][EDGE]  +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][EDGE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][MIDDLE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][MIDDLE]  +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][MIDDLE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] + dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][SIDE_L]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][ISinContigsIds] = dict_info_new[(newStart, newEnd, newISid, newOrient)][ISinContigsIds] + dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][ISinContigsIds]

		del dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)]

			 
def updateMergedSiteInfo(dict_info_new, dict_info_prev, newStart, newEnd, prevStart, prevEnd, oldStart, oldEnd, boolCalcFromContig, directContig, readMapping, edge, middle, side_s, side_l, IStype, orient, IStype_prev, orient_prev):

	print ('New ' + str(newStart) + ':' + str(newEnd) + ':' + IStype + ':' + orient + '\t Prev ' + str(prevStart) + ":" + str(prevEnd) + ":" + IStype_prev + ':' + orient_prev)

	if (newStart, newEnd, IStype, orient) not in dict_info_new: 
		dict_info_new[(newStart, newEnd, IStype, orient)] = {}
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_REF] = boolCalcFromContig
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_CONTIG] = directContig
		dict_info_new[(newStart, newEnd, IStype, orient)][READ_MAPPING] = readMapping
		dict_info_new[(newStart, newEnd, IStype, orient)][EDGE] = edge
		dict_info_new[(newStart, newEnd, IStype, orient)][MIDDLE] = middle
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_L] = side_l
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_S] = side_s

		print ('(1) New ' + str(newStart) + ':' + str(newEnd) + ':' + IStype + ':' + orient + '\t Prev ' + str(prevStart) + ":" + str(prevEnd) + ":" + IStype_prev + ':' + orient_prev)
	else: 
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_REF] = dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_REF] + boolCalcFromContig
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_CONTIG] = dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_CONTIG] + directContig
		dict_info_new[(newStart, newEnd, IStype, orient)][READ_MAPPING] = dict_info_new[(newStart, newEnd, IStype, orient)][READ_MAPPING] + readMapping
		dict_info_new[(newStart, newEnd, IStype, orient)][EDGE] = dict_info_new[(newStart, newEnd, IStype, orient)][EDGE] + edge
		dict_info_new[(newStart, newEnd, IStype, orient)][MIDDLE] = dict_info_new[(newStart, newEnd, IStype, orient)][MIDDLE] + middle
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_L] = dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_L] + side_l
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_S] = dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_S] + side_s

		print ('(2) New ' + str(newStart) + ':' + str(newEnd) + ':' + IStype + ':' + orient + '\t Prev ' + str(prevStart) + ":" + str(prevEnd) + ":" + IStype_prev + ':' + orient_prev)
	
	if prevStart != None and prevEnd != None and not (prevStart == newStart and prevEnd == newEnd and IStype == IStype_prev and orient == orient_prev): 
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_REF] = dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_REF] +  dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][DIRECT_TO_REF]
		dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_CONTIG] = dict_info_new[(newStart, newEnd, IStype, orient)][DIRECT_TO_CONTIG] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][DIRECT_TO_CONTIG]
		dict_info_new[(newStart, newEnd, IStype, orient)][READ_MAPPING] = dict_info_new[(newStart, newEnd, IStype, orient)][READ_MAPPING] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][READ_MAPPING]
		dict_info_new[(newStart, newEnd, IStype, orient)][EDGE] = dict_info_new[(newStart, newEnd, IStype, orient)][EDGE] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][EDGE]
		dict_info_new[(newStart, newEnd, IStype, orient)][MIDDLE] = dict_info_new[(newStart, newEnd, IStype, orient)][MIDDLE] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][MIDDLE]
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_L] = dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_L] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][SIDE_L]
		dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_S] = dict_info_new[(newStart, newEnd, IStype, orient)][SIDE_S] + dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)][SIDE_S] 


		print ('(3) New ' + str(newStart) + ':' + str(newEnd) + ':' + IStype + ':' + orient + '\t Prev ' + str(prevStart) + ":" + str(prevEnd) + ":" + IStype_prev + ':' + orient_prev)
		del dict_info_new[(prevStart, prevEnd, IStype_prev, orient_prev)] 


	"""
	if (newStart, newEnd) not in dict_info_new:
		dict_info_new[(newStart, newEnd)] = {}
		dict_info_new[(newStart, newEnd)][NUM_INF] = dict_info_prev[(prevStart, prevEnd)][NUM_INF]
		dict_info_new[(newStart, newEnd)][COMBINED_SCORE] = dict_info_prev[(prevStart, prevEnd)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd)][SIDE_S] = dict_info_prev[(prevStart, prevEnd)][SIDE_S]
		dict_info_new[(newStart, newEnd)][SIDE_L] = dict_info_prev[(prevStart, prevEnd)][SIDE_L]
        
	else: 
		dict_info_new[(newStart, newEnd)][NUM_INF] = dict_info_new[(newStart, newEnd)][NUM_INF] +  dict_info_prev[(prevStart, prevEnd)][NUM_INF]
		dict_info_new[(newStart, newEnd)][COMBINED_SCORE] = dict_info_new[(newStart, newEnd)][COMBINED_SCORE]  +  dict_info_prev[(prevStart, prevEnd)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd)][SIDE_S] = dict_info_new[(newStart, newEnd)][SIDE_S] +  dict_info_prev[(prevStart, prevEnd)][SIDE_S]
		dict_info_new[(newStart, newEnd)][SIDE_L] = dict_info_new[(newStart, newEnd)][SIDE_L] + dict_info_prev[(prevStart, prevEnd)][SIDE_L]

	if oldStart != None and oldEnd != None and not (newStart == oldStart and newEnd == oldEnd): 
		# print ("Now also add, this info, then del key.")

		dict_info_new[(newStart, newEnd)][NUM_INF] = dict_info_new[(newStart, newEnd)][NUM_INF] +  dict_info_new[(oldStart, oldEnd)][NUM_INF]
		dict_info_new[(newStart, newEnd)][COMBINED_SCORE] = dict_info_new[(newStart, newEnd)][COMBINED_SCORE]  +  dict_info_new[(oldStart, oldEnd)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd)][SIDE_S] = dict_info_new[(newStart, newEnd)][SIDE_S] +  dict_info_new[(oldStart, oldEnd)][SIDE_S]
		dict_info_new[(newStart, newEnd)][SIDE_L] = dict_info_new[(newStart, newEnd)][SIDE_L] + dict_info_new[(oldStart, oldEnd)][SIDE_L]

		del dict_info_new[(oldStart, oldEnd)] 
	"""

def isAnyOverlap(start1, end1, start2, end2, th):

	s1 = start1; e1 = end1; s2 = start2; e2 = end2;

	if start1 > end1:
		s1 = end1
		e1 = start1

	if start2 > end2:
		s2 = end2
		e2 = start2


	if (s1 >= (s2 - th) and s1 <= (e2 + th)) or (e1 >= (s2 - th) and e1 <= (e2 + th)) or (s2 >= (s1 - th) and s2 <= (e1 + th)) or (e2 >= (s1 - th) and e2 <= (e1 + th)):
		return True

	return False


def determineISforRef(dict_contigToRef, dict_IStoContig): 

	list_ISinRef_calc = list() # List of arrLines in gff format

	for refInfo in dict_contigToRef: 
		for (contigId, refStart, refEnd, orient_contigToRef, contigStart, contigEnd, alignLen) in dict_contigToRef[refInfo]: 
			findIn_IStoContig(contigId, contigStart, contigEnd, dict_IStoContig, list_ISinRef_calc, refStart, refEnd, orient_contigToRef, alignLen, refInfo)

	return list_ISinRef_calc


def findIn_IStoContig(contigId_toRef, contigStart_toRef, contigEnd_toRef, dict_IStoContig, list_ISinRef_calc, refStart, refEnd, orient_contigToRef, alignLen_contigToRef, refInfo):
	

	for contigId_toIS in dict_IStoContig: 
		# print(contigId_IS + "\t" + contigId_toRef)
		if (contigId_toRef == contigId_toIS):  
		
			for (contigId_gff, progName, seqFeatType, contigStart_IS, contigEnd_IS, score, orient_IStoContig, phase, attributes) in dict_IStoContig[contigId_toIS]: 

				contigStart_IS = int(contigStart_IS)
				contigEnd_IS = int(contigEnd_IS)     
			
				if ((contigStart_IS >= contigStart_toRef and contigStart_IS <= contigEnd_toRef) or (contigEnd_IS >= contigStart_toRef and contigEnd_IS <= contigEnd_toRef) or (contigStart_toRef >= contigStart_IS and contigStart_toRef <= contigEnd_IS) or (contigEnd_toRef >= contigStart_IS and contigEnd_toRef <= contigEnd_IS)): 



				
					# print(str(refInfo) + "\t" + orient_contigToRef + '\t' + contigId_toRef)
					(calcStart_IStoRef, calcEnd_IStoRef) = convertToRefPos(refStart, refEnd, contigStart_toRef, contigEnd_toRef, contigStart_IS, contigEnd_IS, orient_contigToRef) 

					calcOrient_IStoRef = calcOrientFor_IStoRef(orient_contigToRef, orient_IStoContig)


					attributes = attributes + ';' + 'isCalcFromContig=True'

					tupToAdd = (refInfo, progName, seqFeatType, calcStart_IStoRef, calcEnd_IStoRef, score, calcOrient_IStoRef, phase, attributes) 

					list_ISinRef_calc.append(tupToAdd) 
					# print (contigId_toRef + "\t" + str(contigStart_toRef) + ":" + str(contigEnd_toRef) + "\t" + str(dict_IStoContig[(contigId_toIS)]) + "\t" + orient_contigToRef + ',' + orient_IStoContig + ' -> ' + calcOrient_IStoRef)
					# print (str(tupToAdd))

				# print ('Found!') 

def calcOrientFor_IStoRef(orient_contigToRef, orient_IStoContig):

	if orient_contigToRef == orient_IStoContig: # i.e. '+,+' or '-,-'
		return '+'
	else: 
		return '-'
	 


def convertToRefPos(refStart, refEnd, contigStart_toRef, contigEnd_toRef, contigStart_toIS, contigEnd_toIS, orient_contigToRef):

	# print ('These calculations: ' + str(refStart) + ":" + str(refEnd) + "\t" + str(contigStart_toRef) + ":" + str(contigEnd_toRef) + "\t" + str(contigStart_toIS) + ":" + str(contigEnd_toIS))
	
	calcStart_contig = -1
	calcEnd_contig = -1

	calcStart_ref = -1  
	calcEnd_ref = -1

	if contigStart_toIS >= contigStart_toRef and contigStart_toIS <= contigEnd_toRef: 
		calcStart_contig = contigStart_toIS
	elif contigStart_toIS < contigStart_toRef: 
		calcStart_contig = contigStart_toRef

	if contigEnd_toIS >= contigStart_toRef and contigEnd_toIS <= contigEnd_toRef: 
		calcEnd_contig = contigEnd_toIS
	elif contigEnd_toIS > contigEnd_toRef: 
		calcEnd_contig = contigEnd_toRef

	if (orient_contigToRef == '+'): 
		calcStart_ref = refStart + ( calcStart_contig - contigStart_toRef)
		calcEnd_ref = refStart + ( calcEnd_contig - contigStart_toRef)
	elif (orient_contigToRef == '-'): 
		calcStart_ref = refEnd - ( calcEnd_contig - contigStart_toRef)
		calcEnd_ref = refEnd - ( calcStart_contig - contigStart_toRef)



	# print ('These calculations: ' + str(refStart) + ":" + str(refEnd) + "\t" + str(contigStart_toRef) + ":" + str(contigEnd_toRef) + "\t" + str(contigStart_toIS) + ":" + str(contigEnd_toIS) + '\t' + str(calcStart_ref) + ':' + str(calcEnd_ref))
	# print ("Positions\t" + str(refStart) + ":" + str(refEnd) + "\t" + str(contigStart_toRef) + ":" + str(contigEnd_toRef) + "\t" + str(contigStart_toIS) + ":" + str(contigEnd_toIS) + "\t" + str(calcStart_contig) + ":" + str(calcEnd_contig) + "\t" + str(calcStart_ref) + ":" + str(calcEnd_ref))

	
	return (calcStart_ref, calcEnd_ref)
	

def load_ISinContig_gff(fn_gff_ISannot):

	dict_IStoContig = {} # dict_{(contigId)} => [(arrLine), ...]
	dict_colors = {} # dict_{ISid} => color

	Col_contigId = 0 
	Col_progName = 1
	Col_contigStart = 3 
	Col_contigEnd = 4 

	with open(fn_gff_ISannot, 'r') as fh: 
		for line in fh:
			if not re.match('^#', line):
				line = line.strip() 

				arr = line.split('\t') 

				if arr[Col_progName] == calcInContig.PROG_NAME: 
					# print (line) 

					contigStart = int(arr[Col_contigStart])
					contigEnd = int(arr[Col_contigEnd])

					if arr[Col_contigId] not in dict_IStoContig: 
						dict_IStoContig[arr[Col_contigId]] = []

					dict_IStoContig[arr[Col_contigId]].append(arr)

					extractISidAndColor(arr, dict_colors)

	return (dict_IStoContig, dict_colors)

def load_direct_IStoRef(fn_directIStoRef, isMerged, isIgnoreOrient, isIgnoreIStype, separator):

	dict_direct = calcInContig.loadDirectIStoContig(fn_directIStoRef, [], 0, 0)

	dict_directMerged = dict() # dict_{(refId, refLen)} => [(ISid, alignPosInRef_start, alignPosInRef_end, orientIStoContig, isContigFlipped)]

	for (refId, refLen) in dict_direct: 

		if (refId, refLen) not in dict_directMerged: 
			dict_directMerged[(refId, refLen)] = [] 

		for (ISid_1, refStart_1, refEnd_1, orient_1, isContigFlipped_1) in dict_direct[(refId, refLen)]: 

			if len(dict_directMerged[(refId, refLen)]) == 0: 
				dict_directMerged[(refId, refLen)].append((ISid_1, refStart_1, refEnd_1, orient_1, isContigFlipped_1))

			else: 
				# Either merge or add new! 
				isUpdated = checkIfMerge_directRefIS(dict_directMerged[(refId, refLen)], isMerged, isIgnoreOrient, isIgnoreIStype, separator, ISid_1, refStart_1, refEnd_1, orient_1, isContigFlipped_1)

				if isUpdated == False: 
					dict_directMerged[(refId, refLen)].append((ISid_1, refStart_1, refEnd_1, orient_1, isContigFlipped_1))
				
	return dict_directMerged 

def identifyISonlyContigs(fn_IStoContig, th_alignDiff_IStoContig):
	# Identify contigs which are completely IS (for exclusion from further analysis)

	list_contigsWithComplIS = [] 

	with open(fn_IStoContig, 'r') as fh:  
		for line in fh: 
			if not re.match('^#', line):
				line = line.strip() 

				arr = line.split('\t') 

				# print (line)
				# print (arr[Col_alignLen] + '\t' + arr[Col_sSeqLen])
				alignDiff = abs(int(arr[Col_alignLen]) - int(arr[Col_sSeqLen]))
				if alignDiff <= th_alignDiff_IStoContig:  
					# print (line)
					list_contigsWithComplIS.append(arr[Col_sSeqId]) 

	print ('List with complete IS')
	print (list_contigsWithComplIS) 
	return (list_contigsWithComplIS)

def getNewISid_directRefIS(ISid_orig, ISid_merged, separator):

	ISid_new = '' 

	arr = ISid_merged.split(separator)
	
	if ISid_orig in arr: 
		ISid_new = ISid_merged 

	else: 
		arr.append(ISid_orig)
		ISid_new = separator.join(sorted(arr))
		# ISid_new = ISid_merged + separator + ISid_orig

	return ISid_new
	

def getISnotMappedToRef(list_allContigs, dict_contigToRef_alignLen, dict_IStoContig, list_contigsWithComplIS):

	print('**************************************************** ISinContig:' + str(len(dict_IStoContig)) + ' ContigToRef:' + str(len(dict_contigToRef_alignLen)))

	# for aContigId in dict_IStoContig: 
		# if aContigId not in list_contigsWithComplIS and (aContigId not in dict_contigToRef_alignLen or isAlignLenDiffHigh(dict_contigToRef_alignLen[aContigId])) : 
		#	print ('Print all IS in this contig (if not covered in alignment). ' + aContigId)
		#	print (dict_IStoContig[aContigId])



		# else: 
			# Check aligned len
		# 	print ('Aligned len ' + str(dict_contigToRef_alignLen[aContigId]))
"""
	for aContigId in list_allContigs: 
		if aContigId not in dict_contigToRef_alignLen: 
			print('ContigId not alignedToRef: ' + aContigId)
			if aContigId in dict_IStoContig:
				print ('IS not aligned to ref:')
				print (dict_IStoContig[aContigId])


"""
def isAlignLenDiffHigh(dict_contigToRef_alignInfo):

	if dict_contigToRef_alignInfo['total'] - dict_contigToRef_alignInfo['covered'] > 1: 
		return True
	
	return False

def load_contigToRef(fn_blastRes_contigToRef, list_contigsWithComplIS, dict_direct, th_maxAlignLenDiff, dict_IStoContig, th_contigToRef_ISonlyLenDiff): 
	## MAY NEED TO FILTER THE ALIGNMENTS. 

	dict_contigsAligned = {} # {contigId} => {total: , covered: , positions: []}

	dict_contigToRef = {} # dict_{(refId, refLen)} => [(ContigId, refStart, refEnd, orient_contigToRef, contigStart, contigEnd, alignLen)]
	
	list_allContigs = [] # [contigs1, contigs2, contigs3, ..., contigsN]

	with open(fn_blastRes_contigToRef, 'r') as fh: 

		contigId_counter = 0
		contigId_prev = None

		for line in fh: 
			line = line.strip() 
			
			if re.match('^# Query: ', line):
				
				contigId = re.sub('^# Query: ', '', line)
				# print (contigId)

				if contigId not in list_allContigs: 
					
					list_allContigs.append(contigId)


			elif not re.match('^#', line):
				
				arr = line.split('\t')

				# print ('ContigId ...? ' + arr[Col_qSeqId])

				if (arr[Col_qSeqLen] not in list_contigsWithComplIS): 

					if arr[Col_qSeqId] != contigId_prev: 
						# print ('New contig encountered ' + arr[Col_qSeqId])
						contigId_counter = 0; 
						contigId_prev = arr[Col_qSeqId]

					refStart = int(arr[Col_sStart])
					refEnd = int(arr[Col_sEnd])

					contigStart = int(arr[Col_qStart])
					contigEnd = int(arr[Col_qEnd])

					orient = '' 
					if (refStart > refEnd and contigStart > contigEnd) or (refStart < refEnd and contigStart < contigEnd): # same direction
						orient = '+'
						# print('\t+')
					else:
						orient = '-'
						# print('\t-')


					if refStart > refEnd:
						tmp = refStart
						refStart = refEnd
						refEnd = tmp

					if contigStart > contigEnd: 
						tmp = contigStart
						contigStart = contigEnd
						contigEnd = tmp 


					## Do some filtering...
					# Filtering if alignment region is exactly the same as ISinRef. 
					if isISonlyAlignedRegion(arr[Col_qSeqId], contigStart, contigEnd, dict_IStoContig, th_contigToRef_ISonlyLenDiff) == True: 
						continue 

					# Filtering: 
					if (arr[Col_qSeqId] not in dict_contigsAligned or not (dict_contigsAligned[arr[Col_qSeqId]]['total'] - dict_contigsAligned[arr[Col_qSeqId]]['covered'] <= th_maxAlignLenDiff)) :  

						contigId_counter = contigId_counter + 1 

						## Add to dict 
						if (arr[Col_sSeqId], arr[Col_sSeqLen]) not in dict_contigToRef: 
							dict_contigToRef[(arr[Col_sSeqId], arr[Col_sSeqLen])] = [] 

						dict_contigToRef[(arr[Col_sSeqId], arr[Col_sSeqLen])].append((arr[Col_qSeqId], refStart, refEnd, orient, contigStart, contigEnd, int(arr[Col_alignLen])))

						# if int(arr[Col_alignLen]) == int(arr[Col_qSeqLen]):

						if arr[Col_qSeqId] not in dict_contigsAligned: 
							dict_contigsAligned[arr[Col_qSeqId]] = {} 
							dict_contigsAligned[arr[Col_qSeqId]]['total'] = int(arr[Col_qSeqLen])
							dict_contigsAligned[arr[Col_qSeqId]]['covered'] = contigEnd - contigStart + 1 
							dict_contigsAligned[arr[Col_qSeqId]]['positions'] = [contigStart, contigEnd]
						else: 
							dict_contigsAligned[arr[Col_qSeqId]]['covered'] = updateContigsPositions(dict_contigsAligned[arr[Col_qSeqId]]['positions'], contigStart, contigEnd) 


					# list_completeAlignment.append(arr[Col_qSeqId])     
					# print ('ADDING THIS THING')

					
						# print ('False') 
							# print (str(refStart) + ":" + str(refEnd) + "\t" + str(contigStart) + ":" + str(contigEnd) + "\t" + orient)

					#else: 
					#    print ('True')
	print ('################')
	print ('Num. of all contigs in contigToRef blastRes ' + str(len(list_allContigs)))
	print ('Num. of dict_contigsAligned ' + str(len(dict_contigsAligned))) 

	
	for refId in dict_contigToRef: 
		print ('Num. of contigsToRef ' + str(len(dict_contigToRef[refId])))

	print ('################')

	for key in dict_contigsAligned: 
		print ('ContigsAligned: ' + str(key) + "\t" + str(dict_contigsAligned[key])) 
		

	return (dict_contigToRef, list_allContigs, dict_contigsAligned)



def updateContigsPositions(list_s_e, contigStart, contigEnd): 

	if contigStart < list_s_e[0]: 
		list_s_e[0] = contigStart

	if contigEnd > list_s_e[1]: 
		list_s_e[1] = contigEnd 

	return (list_s_e[1] - list_s_e[0] + 1)


def isAlignmentInDirRefIS(newRefStart, newRefEnd, dict_direct): # Is checking for overlaps now ... 
	"""
	for refInfo in dict_direct: 
		for (contigId, refStart, refEnd, orient, isContigFlipped) in  dict_direct[refInfo]: 
			
			
			if (newRefStart >= refStart and newRefStart <= refEnd) or (newRefEnd >= refStart and newRefEnd <= refEnd): 
				# print (str(newRefStart) + ":" + str(newRefEnd) + "\t" + str(refStart) + ":" + str(refEnd))
				return True
	"""

	return False      

def calcIfInDirect(dict_direct, IStoRef_refInfo, IStoRef_start, IStoRef_end, IStoRef_orient, IStoRef_attributes): 
	print ("# IStoRef")
	print (IStoRef_refInfo)
	print (IStoRef_attributes)


def addCounts_ISidsInContigs(dict_, list_uniqISids): 

	for uniqId in list_uniqISids: 
		if uniqId not in dict_: 
			dict_[uniqId] = 0

		dict_[uniqId] = dict_[uniqId] + 1

def calcAndPrintISnotMappedToRef(dict_ISidsInContigs_counts, dict_IStoContig):

	### TO DELETE!!!! 
	del dict_ISidsInContigs_counts['1']
	del dict_ISidsInContigs_counts['2']
	del dict_ISidsInContigs_counts['3']
	del dict_ISidsInContigs_counts['4']
	del dict_ISidsInContigs_counts['5']
	del dict_ISidsInContigs_counts['6']
	del dict_ISidsInContigs_counts['7']
	del dict_ISidsInContigs_counts['8']
	del dict_ISidsInContigs_counts['9']
	del dict_ISidsInContigs_counts['10']

	list_allISids = [] 

	for contigId in dict_IStoContig: 
		for arrLine in dict_IStoContig[contigId]: 
			uniqId = extractUniqId(arrLine[Col_gff3_info])
			
			if uniqId not in list_allISids: 
				list_allISids.append(uniqId)
			
			

	print('================================= ' + str(len(list_allISids)) + ' ' + str(len(dict_ISidsInContigs_counts)))

	list_uniqISids_notInRef = list(set(list_allISids) - set(dict_ISidsInContigs_counts))


	print ('Diff ids are: ')
	print (list_uniqISids_notInRef)
	
	for uniqISid in OrderedDict(sorted(dict_ISidsInContigs_counts.items(), key=lambda t: int(t[0]))): 
		print (str(uniqISid) + '\t' + str(dict_ISidsInContigs_counts[uniqISid]))

def extractUniqId(infoStr): 
	arr = infoStr.split(';')

	for val in arr: 
		if re.match('^uniqId\=', val): 
			uniqId = re.sub('^uniqId\=', '', val)

	return uniqId

#################################### MAIN 
def main(): 

	sys.stderr.write("\nNotice: If you are ussing the --isMerged setting, please run mergeLocalCounts.py setting isIgnoreOrient and isIgnoreIStype to False. Otherwise only --isMerged counts (and isIgnoreOrient counts) for the ref will be overestimated here.\n\n")

	parser = argparse.ArgumentParser(description='Calculate IS insertion position, ISid, and orientation in comparison to reference.')


	parser.add_argument('--IStoContig', nargs=1, required=True, help='Blast results of IS to contigs.')
	parser.add_argument('--th_alignDiff_IStoContig', nargs=1, default=[5], type=int, help="Maximum alignment difference between align-length and contig-length to designate a contig as completeIS. This contig's alignment to the reference is then excluded. Default=5.")


	parser.add_argument('--directIStoRef', nargs=1, required=True, help='Blast results of IS to reference.')

	parser.add_argument('--contigToRef', nargs=1, required=True, help='Blast results of contig to reference.')


	parser.add_argument('--IS_annotation', nargs=1, required=True, help='IS annotations file in gff3 format.')

	parser.add_argument('--fn_out', nargs=1, required=True, help='IS annotations for reference in gff3 format.')
	## parser.add_argument('--fn_out_notInRef', nargs=1, required=True, help='IS annotations not in reference in gff3 format.')

	parser.add_argument('--isMerged', nargs=1, required=False, help='Only print the mapped IS from contigs (from the merged set). Default is False.', default=['False'])
	parser.add_argument('--th_merge', nargs=1, required=False, help='Threshold (i.e. sequence length in base pairs) to merge overlapping IS. Value only used in conjection with --isMerged. Default=20.', default=[20])
	parser.add_argument('--ignoreOrient', nargs=1, required=False, help='Values as True or False - to igore or-else check orientation, repsectively. Ignored orientation is appended as orient1:orient2. Default=False.', default=["False"])
	parser.add_argument('--ignoreIStype', nargs=1, required=False, help='Values as True or False - to igore or-else check IStype, repsectively. Ignored IStypes are appended as IStype1:IStype2:..:IStypeN. This automatically sets isIgnoreOrient to "True". Default=False.', default=["False"])
	parser.add_argument('--separator', required=False, default=[':'], help='Separate unresolvable IS and orientations. Default=":". Choose a separator that is not an alphabet in the IS identifier.')

	parser.add_argument('--th_maxContigToRef_AlignLenDiff', required=False, default=[10], help='Total len (of contigs) - AlignedLen (of contigs) <= this_threshold. Default = 10.')

	parser.add_argument('--th_contigToRef_ISonlyLenDiff', required=False, default=[50], help='Alignment length difference between contigToRef and IStoContig - to remove alignments of just the IS in contigs aligning to the reference contigs.')

	parser.add_argument('--th_directIS_overlap', nargs=1, required=False, help='Threshold to check if identified IS overlaps with those found in the directIS. Default=20.', default=[20])
	args = parser.parse_args()

	calcPosToRef(args.contigToRef[0], args.IStoContig[0], args.th_alignDiff_IStoContig[0], args.IS_annotation[0], args.directIStoRef[0], args.fn_out[0], args.isMerged[0], int(args.th_merge[0]), args.ignoreOrient[0], args.ignoreIStype[0], args.separator[0], int(args.th_directIS_overlap[0]), int(args.th_maxContigToRef_AlignLenDiff[0]), int(args.th_contigToRef_ISonlyLenDiff[0]))

if __name__ == '__main__':
	main()