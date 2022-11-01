#!/usr/bin/python3

import argparse
import re
import sys
import calcInContig_posISOrient



EDGE = 'Edge'
MIDDLE = 'Middle'
ESTIMATE = 'Total_estimate'

SITES = 'sites'
INFO = 'information'
NUM_INF = 'numInf'
COMBINED_SCORE = 'combinedScore'
SIDE_S = 'side_s'
SIDE_L = 'side_l'

GFF3_VERSION = '##gff-version 3.1.26'


############################################ TOP_LVL
def getTheFinalCountOfIS(fn_ISinGff, th_forMergingOverlaps, isIgnoreOrient, th_toCountAsEdge, isIgnoreIStype, fnOut_estimates, separator, fnOut_gff3_merged, isolateId, fnOut_estimates_singleRow):

	(dict_allData, dict_seqReqInfo, dict_colors) = loadTheData(fn_ISinGff)

	dict_reducedSet1 = reduceSet_1(dict_allData, dict_seqReqInfo) # Merges an overlapping (without using the threshold), no final integrity check!

	dict_reducedSet2 = reduceSet_2(dict_reducedSet1, th_forMergingOverlaps) # Merge adjacent (minimally), and using the specified threshold
	
	## Do an integrity check! (EXIT if check fails)
	#checkIntegrity(dict_reducedSet2, th_forMergingOverlaps)


	if isIgnoreIStype == "True":
		# print("Calculating reduced set ")
		dict_reducedSet3 = reduceSet_3(dict_reducedSet2, th_forMergingOverlaps, separator)
		dict_reducedSet4 = reduceSet_4(dict_reducedSet3, th_forMergingOverlaps, separator)


		(dict_estimates, dict_edgeOrMiddleAnnot) = doAnEdgeInnerCount(dict_reducedSet4, dict_seqReqInfo, th_toCountAsEdge)
		printEstimateToFile(dict_estimates, fnOut_estimates)
		printEstimateToFile_singleRow(dict_estimates, fnOut_estimates_singleRow, isolateId)
		printGFF3_merged(dict_reducedSet4, separator, dict_edgeOrMiddleAnnot, fnOut_gff3_merged, dict_seqReqInfo, dict_colors)

	elif isIgnoreOrient == "True":
		# print ('Calc a reducedSet3')
		dict_reducedSet3 = reduceSet_3(dict_reducedSet2, th_forMergingOverlaps, separator)
		(dict_estimates, dict_edgeOrMiddleAnnot) = doAnEdgeInnerCount(dict_reducedSet3, dict_seqReqInfo, th_toCountAsEdge)
		printEstimateToFile(dict_estimates, fnOut_estimates)
		printEstimateToFile_singleRow(dict_estimates, fnOut_estimates_singleRow, isolateId)
		printGFF3_merged(dict_reducedSet3, separator, dict_edgeOrMiddleAnnot, fnOut_gff3_merged, dict_seqReqInfo, dict_colors)
	
	else:
		(dict_estimates, dict_edgeOrMiddleAnnot) = doAnEdgeInnerCount(dict_reducedSet2, dict_seqReqInfo, th_toCountAsEdge)
		printEstimateToFile(dict_estimates, fnOut_estimates)
		printEstimateToFile_singleRow(dict_estimates, fnOut_estimates_singleRow, isolateId)
		printGFF3_merged(dict_reducedSet2, separator, dict_edgeOrMiddleAnnot, fnOut_gff3_merged, dict_seqReqInfo, dict_colors)

	


	"""
	for key in dict_allData:
		print ("# " + str(key))
		print (dict_allData[key])
		print ()
	"""

	"""
	for contigId in dict_reducedSet1:
		print (contigId)

		for key in dict_reducedSet1[contigId]:
			print ('\t' + str(key) + str(dict_reducedSet1[contigId][key]))

		print ()
	"""

	"""
	for key in dict_reducedSet2:
		print (key)

		for aSite in dict_reducedSet2[key]:
			print ('\t' + str(aSite) + '\t' + str(dict_reducedSet2[key][aSite]))
	"""

	"""
	for key in dict_reducedSet3:
		print (key)

		for aSite in dict_reducedSet3[key]:
			print ('\t' + str(aSite) + '\t' + str(dict_reducedSet3[key][aSite]))
	"""

	"""
	for key in dict_reducedSet4:
		print (key)

		for aSite in dict_reducedSet4[key]:
			print ('\t' + str(aSite) + '\t' + str(dict_reducedSet4[key][aSite]))
	"""
			
	"""
	for contigId in dict_reducedSet4:


		for (start, end, ISid, orient) in dict_reducedSet4[contigId][]:
			print (contigId + '\t' + 'WiIS' + '\t' + 'insertion_sequence' + '\t' + str(start) + '\t' + str(end) + '\t' + 'score' + '\t' + orient + '\t' + '.' + '\t' + 'attributes')
	"""

############################################ AUX - OUTPUT
def printGFF3_merged(dict_reducedSet4, separator, dict_edgeOrMiddleAnnot, fnOut_gff3_merged, dict_seqReqInfo, dict_colors):

	fh_out = open(fnOut_gff3_merged, 'w+')
	# Print the file header
	fh_out.write(GFF3_VERSION + '\n') 
	for contigId in dict_reducedSet4:

		fh_out.write('##sequence-region ' + contigId + ' ' + str(dict_seqReqInfo[contigId][0]) + ' ' +  str(dict_seqReqInfo[contigId][1]) + '\n') 
		pass

	# Print the data

	uniqId = 0
	for contigId in dict_reducedSet4:
		for (start, end, ISid, orient) in dict_reducedSet4[contigId][SITES]:

			if ISid not in dict_colors: 
				dict_colors[ISid] = calcInContig_posISOrient.generateRandomColor() 

			fh_out.write (contigId + '\t' + 'WiIS' + '\t' + 'insertion_sequence' + '\t' + str(start) + '\t' + str(end) + '\t' + '.')
			if separator not in orient:
				fh_out.write('\t' + orient)
			else:
				fh_out.write('\t' + '.')

			fh_out.write('\t' + '.' + '\t' + 'Name=' + ISid)

			dict_info = dict_reducedSet4[contigId][INFO][(start, end, ISid, orient)] 
			
			fh_out.write(';num_Inf=' + str(dict_info[NUM_INF]))
			fh_out.write(';combined_Depth=' + str(dict_info[COMBINED_SCORE]))
			fh_out.write(';side_s=' + str(dict_info[SIDE_S]))
			fh_out.write(';side_l=' + str(dict_info[SIDE_L]))


			fh_out.write(';positionInContig=' + ','.join(dict_edgeOrMiddleAnnot[contigId][start, end, ISid, orient]))

			fh_out.write(';color=' + dict_colors[ISid])

			fh_out.write(';uniqId=' + str(uniqId))
			uniqId = uniqId + 1

			fh_out.write('\n')


	fh_out.close() 

def printEstimateToFile_singleRow(dict_eachIScount, fnOut_estimates, isolateId): 

	list_IStypes = sorted(dict_eachIScount)


	fh_out = open(fnOut_estimates, 'w+')

	# Write the header
	fh_out.write('\t')
	fh_out.write('#Assembly') 

	for IStype in list_IStypes:
		for calculatedType in [EDGE, MIDDLE, ESTIMATE]:
			fh_out.write('\t' + IStype + "_WaIS:" + calculatedType);

	for calculatedType in [EDGE, MIDDLE, ESTIMATE]:
		fh_out.write('\t' + "Total_WaIS:" + calculatedType);

	# fh_out.write('\t')
	# fh_out.write('Total')
	fh_out.write('\n')
	
	# The output content 
	fh_out.write(isolateId)
	fh_out.write('\t')
	for IStype in list_IStypes:
		total_posInContig = 0
		for calculatedType in [EDGE, MIDDLE]:
			fh_out.write('\t')
			if calculatedType in dict_eachIScount[IStype]:
				fh_out.write(str(dict_eachIScount[IStype][calculatedType]))
				if calculatedType == EDGE: 
					total_posInContig = total_posInContig + (dict_eachIScount[IStype][calculatedType]/2)
				elif calculatedType == MIDDLE: 
					total_posInContig = total_posInContig + dict_eachIScount[IStype][calculatedType]
			else:
				fh_out.write('0')

		fh_out.write('\t' + str(total_posInContig))

	totalEdge = 0; totalMiddle = 0; totalEst = 0; 
	for IStype in list_IStypes:

		if EDGE in dict_eachIScount[IStype]:
			totalEdge = totalEdge + dict_eachIScount[IStype][EDGE]

		if MIDDLE in dict_eachIScount[IStype]: 
			totalMiddle = totalMiddle + dict_eachIScount[IStype][MIDDLE]

	totalEst = totalMiddle + (totalEdge/2)

	fh_out.write('\t' + str(totalEdge))
	fh_out.write('\t' + str(totalMiddle))
	fh_out.write('\t' + str(totalEst))

	fh_out.write('\n')

	fh_out.close() 

def printEstimateToFile(dict_eachIScount, fnOut_estimates):

	fh_out = open(fnOut_estimates, 'w+')

	fh_out.write('#################### Assemby' + "\n")
	list_IStypes = sorted(dict_eachIScount)

	fh_out.write('Tool')
	for IStype in list_IStypes:
		fh_out.write('\t' + IStype);
	fh_out.write('\tTotal\n')


	for posInContig in [EDGE, MIDDLE]:
		fh_out.write('WiIS:' + posInContig)

		total_posInContig = 0

		for IStype in list_IStypes:

			fh_out.write("\t")
			if posInContig in dict_eachIScount[IStype]:
				fh_out.write(str(dict_eachIScount[IStype][posInContig]))
				total_posInContig = total_posInContig + dict_eachIScount[IStype][posInContig]
			else:
				fh_out.write('0')

			# print (IStype + '\t' + posInContig + "\t" + str(dict_eachIScount[IStype][posInContig]))

		fh_out.write("\t" + str(total_posInContig) + '\n')


	fh_out.write('WiIS:' + ESTIMATE)
	totalEst = 0
	for IStype in list_IStypes:

		edge = 0 if EDGE not in dict_eachIScount[IStype] else dict_eachIScount[IStype][EDGE]
		middle = 0 if MIDDLE not in dict_eachIScount[IStype] else dict_eachIScount[IStype][MIDDLE]

		estimate = (edge/2) + middle
		fh_out.write('\t' + str(estimate))

		totalEst = totalEst + estimate

	fh_out.write('\t' + str(totalEst) + '\n')

	fh_out.close()


############################################ AUX - REDUCED SETS
def doAnEdgeInnerCount(dict_reducedSet3, dict_seqReqInfo, th_toCountAsEdge):

	dict_eachIScount = dict() # dict_
	dict_edgeOrMiddleAnnot = dict() # dict_[contigId] => (start, end) => 'Edge|Middle'

	countEdge = 0; countMiddle = 0

	for contigId in dict_reducedSet3:

		for (start, end, ISid, orient) in dict_reducedSet3[contigId][SITES]:

				(contigStart, contigEnd) = dict_seqReqInfo[contigId]


				contigStart_withTh = (contigStart + th_toCountAsEdge) if (contigStart + th_toCountAsEdge) < contigEnd else contigEnd
				contigEnd_withTh = (contigEnd - th_toCountAsEdge) if (contigEnd - th_toCountAsEdge) > contigStart else contigStart


				if isAnyOverlap(start, end, contigStart, contigStart_withTh, 0) or isAnyOverlap(start, end, contigEnd_withTh, contigEnd, 0):

					# print(contigId + '\t' +  EDGE + ": " + str((start, end, ISid, orient)))
					countEdge = countEdge + 1

					addToEachISCount(dict_eachIScount, ISid, EDGE )

					if contigId not in dict_edgeOrMiddleAnnot: 
						dict_edgeOrMiddleAnnot[contigId] = {}
					if (start, end, ISid, orient) not in dict_edgeOrMiddleAnnot[contigId]: 
						dict_edgeOrMiddleAnnot[contigId][(start, end, ISid, orient)] = []
					dict_edgeOrMiddleAnnot[contigId][(start, end, ISid, orient)].append(EDGE)


				else:

					# print(contigId + "\t" + MIDDLE + ": " + str((start, end, ISid, orient)))
					countMiddle = countMiddle + 1

					addToEachISCount(dict_eachIScount, ISid, MIDDLE)
					
					if contigId not in dict_edgeOrMiddleAnnot: 
						dict_edgeOrMiddleAnnot[contigId] = {}
					if (start, end, ISid, orient) not in dict_edgeOrMiddleAnnot[contigId]: 
						dict_edgeOrMiddleAnnot[contigId][(start, end, ISid, orient)] = []
					dict_edgeOrMiddleAnnot[contigId][(start, end, ISid, orient)].append(MIDDLE)	


	return (dict_eachIScount, dict_edgeOrMiddleAnnot)

def addToEachISCount(dict_eachIScount, ISid, posInContig):

	if ISid not in dict_eachIScount:
		dict_eachIScount[ISid] = {}

	if posInContig not in dict_eachIScount[ISid]:
		dict_eachIScount[ISid][posInContig] = 0

	dict_eachIScount[ISid][posInContig] = dict_eachIScount[ISid][posInContig] + 1



def checkIntegrity(dict_reducedSet2, th):
	START = 0; END = 1; IS = 2; ORIENT = 3

	for refId in dict_reducedSet2:
		list_sites = dict_reducedSet2[refId]

		# print("Loop: " + str(len(list_sites)))
		for i in range(0, len(list_sites)):
			for j in range (i+1, len(list_sites)):
				if i != j:

					if list_sites[i][IS] == list_sites[j][IS] and list_sites[i][ORIENT] == list_sites[j][ORIENT] and isAnyOverlap(list_sites[i][START], list_sites[i][END], list_sites[j][START], list_sites[j][END], th):
						sys.stderr.write("Error: overlap detected in two sites " + str(list_sites[i][START]) + ':' + str(list_sites[i][END]) + '\t' + str(list_sites[j][START]) + ':' + str(list_sites[j][END]) + '\t' + list_sites[j][IS] + '\t' + list_sites[j][ORIENT])
						sys.exit()

					# print (str(i) + ' ' + str(j))

	# pass



def reduceSet_4(dict_reducedSet1, th_forMergingOverlaps, separator):

	dict_reducedSet2 = dict() # dict_reducedSet2[refId|contigId] = [(startPos|min, endPos|max, ISid, 'orient1:orient2')]

	for key in dict_reducedSet1:
		for (start_set1, end_set1, ISid_set1, orient_set1) in dict_reducedSet1[key][SITES]:

			if key not in dict_reducedSet2:
				dict_reducedSet2[key] =  {SITES: [], INFO: {}}
				dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
				addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			else: # Loop-through, check - merge or create new

				isOverlapFound = False
				for idx, (start_set2, end_set2, ISid_set2, orient_set2) in enumerate(dict_reducedSet2[key][SITES]):
					if isAnyOverlap(start_set1, end_set1, start_set2, end_set2, th_forMergingOverlaps):

						# dict_redu
						# print ('Do the merge!')

						newStart = start_set2 if (start_set2 <= start_set1) else start_set1
						newEnd = end_set2 if end_set2 >= end_set1 else end_set1


						orient = separator.join(list((set(orient_set1.split(separator) + orient_set2.split(separator)))))


						""" 
						ISid = ISid_set2

						# print (ISid_set1 + '\t' + ISid_set2)
						if separator in ISid_set2:
							arr_ = ISid_set2.split(separator)
							# print (arr_)
							if ISid_set1 not in arr_:
								arr_.append(ISid_set1)
								ISid = separator.join(sorted(arr_))
								# ISid = ISid + separator + ISid_set1
						elif ISid_set1 != ISid_set2:
							ISid = separator.join(sorted([ISid_set2,  ISid_set1]))
							# ISid = ISid_set2 + separator + ISid_set1
						"""

						ISid = separator.join(list((set(ISid_set1.split(separator) + ISid_set2.split(separator)))))	

						dict_reducedSet2[key][SITES][idx] = (newStart, newEnd, ISid, orient)
						addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2, ISid, orient, ISid_set1, orient_set1, ISid_set2, orient_set2)

						isOverlapFound = True
						break

				if isOverlapFound == False:
					dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
					addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)


			# 	if isAnyOverlapß
			# print ('\t' + str())

# 	pass

	return dict_reducedSet2

def reduceSet_3(dict_reducedSet1, th_forMergingOverlaps, separator):

	dict_reducedSet2 = dict() # dict_reducedSet2[refId|contigId] = {'sites' => [(startPos|min, endPos|max, ISid, 'orient1:orient2')]; 'information' => {numInf: '0', combinedScore: '', side_s: '', side_l: ''}]

	for key in dict_reducedSet1:
		for (start_set1, end_set1, ISid_set1, orient_set1) in dict_reducedSet1[key][SITES]:

			if key not in dict_reducedSet2:
				dict_reducedSet2[key] = {SITES: [], INFO: {}}
				dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))

				addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			else: # Loop-through, check - merge or create new

				isOverlapFound = False
				for idx, (start_set2, end_set2, ISid_set2, orient_set2) in enumerate(dict_reducedSet2[key][SITES]):
					if ISid_set1 == ISid_set2 and isAnyOverlap(start_set1, end_set1, start_set2, end_set2, th_forMergingOverlaps):

						# dict_redu
						# print ('Do the merge!')

						newStart = start_set2 if (start_set2 <= start_set1) else start_set1
						newEnd = end_set2 if end_set2 >= end_set1 else end_set1


						orient = separator.join(list((set(orient_set1.split(separator) + orient_set2.split(separator)))))

						dict_reducedSet2[key][SITES][idx] = (newStart, newEnd, ISid_set2, orient)
						addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2, ISid_set1, orient, ISid_set1, orient_set1, ISid_set2, orient_set2)

						isOverlapFound = True
						break

				if isOverlapFound == False:
					dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
					addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)


			# 	if isAnyOverlapß
			# print ('\t' + str())

# 	pass

	return dict_reducedSet2






def reduceSet_2(dict_reducedSet1, th_forMergingOverlaps):

	dict_reducedSet2 = dict() # # dict_reducedSet2[refId|contigId] = {'sites' => [(startPos|min, endPos|max, ISid, 'orient1:orient2')]; 'information' => {(startMin, endMax)} => {numInf: '0', combinedScore: '', side_s: '', side_l: ''}]


	for key in dict_reducedSet1:
		for (start_set1, end_set1, ISid_set1, orient_set1) in dict_reducedSet1[key][SITES]:

			if key not in dict_reducedSet2:
				dict_reducedSet2[key] = {SITES: [], INFO: {}}
				dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))

				addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

				# print ("Here! " + str(dict_reducedSet1[key][INFO][(start_set1, end_set1)]))
				# info = dict_reducedSet1[key][INFO][(start_set1, end_set1)]
				# addSupportingInfo(dict_reducedSet2, key, start_set1, end_set1, score, side_s, side_l, None, None)

			else: # Loop-through, check - merge or create new

				isOverlapFound = False
				for idx, (start_set2, end_set2, ISid_set2, orient_set2) in enumerate(dict_reducedSet2[key][SITES]):
					if ISid_set1 == ISid_set2 and orient_set1 == orient_set2 and isAnyOverlap(start_set1, end_set1, start_set2, end_set2, th_forMergingOverlaps):

						# dict_redu
						# print ('Do the merge!')

						newStart = start_set2 if (start_set2 <= start_set1) else start_set1
						newEnd = end_set2 if end_set2 >= end_set1 else end_set1


						dict_reducedSet2[key][SITES][idx] = (newStart, newEnd, ISid_set2, orient_set2)
						addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], newStart, newEnd, start_set1, end_set1, start_set2, end_set2, ISid_set1, orient_set1, ISid_set1, orient_set1, ISid_set2, orient_set2)

						isOverlapFound = True
						break

				if isOverlapFound == False:
					dict_reducedSet2[key][SITES].append((start_set1, end_set1, ISid_set1, orient_set1))
					addSupportingInfo_reduced(dict_reducedSet2[key][INFO], dict_reducedSet1[key][INFO], start_set1, end_set1, start_set1, end_set1, None, None, ISid_set1, orient_set1, ISid_set1, orient_set1, None, None)

			# 	if isAnyOverlap
			# print ('\t' + str())

# 	pass

	return dict_reducedSet2

def reduceSet_1(dict_allData, dict_seqReqInfo):
	dict_reducedSet = dict() # dict_reducedSet2[refId|contigId] = {'sites' => [(startPos|min, endPos|max, ISid, 'orient1:orient2')]; 'information' => {(startMin, endMax)} => {numInf: '0', combinedScore: '', side_s: '', side_l: ''}]


	list_allISsites = sorted(dict_allData.keys(), key=lambda x: x[1])

	for (refId, startPos, endPos, ISid, orient, isCalcFromContig) in list_allISsites:




		list_info = dict_allData[(refId, startPos, endPos, ISid, orient, isCalcFromContig)]

		for (modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, score) in list_info:



			# print ('score: ' + str(score) + ' side_s: ' + str(side_s) + ' side_l:' + str(side_l) )

			checkAndAddASite(dict_reducedSet, refId, score, startPos, endPos, minPos, maxPos, ISid, orient, side_s, side_l)

			# if (score != 'Inf'):

			# 	print (str(score) + ' ' + str(startPos) + ':' + str(endPos) + " " + str(minPos) + ":" + str(maxPos))

		"""
		if refId not in dict_reducedSet:
			dict_reducedSet[refId] = []
			dict_reducedSet[refId].append((startPos, endPos, ISid, orient))

		else:
			# Check, then add

		"""

	return dict_reducedSet

def checkAndAddASite(dict_reducedSet, refId, score, startPos, endPos, minPos, maxPos, ISid, orient, side_s, side_l):

	# Input data handling
	startPos = int(startPos)
	endPos = int(endPos)

	if startPos > endPos:
		tmp = startPos
		startPos = endPos
		endPos = tmp

	if score != 'Inf':
		minPos = int(minPos)
		maxPos = int(maxPos)

		tmp = minPos
		minPos = maxPos
		maxPos = tmp

	if score == 'Inf':
		start_min = startPos
		end_max = endPos
	else:
		start_min = startPos if startPos <= minPos else minPos
		end_max = endPos if endPos >= maxPos else maxPos

	# print ('After: ' + str(start_min) + ':' + str(end_max))

	if refId not in dict_reducedSet:
		dict_reducedSet[refId] = {SITES: [], INFO: {}} #


	if len(dict_reducedSet[refId]) == 0: # Adding the first site
		dict_reducedSet[refId][SITES].append((start_min, end_max, ISid, orient))
		addSupportingInfo(dict_reducedSet, refId, start_min, end_max, score, side_s, side_l, None, None, ISid, orient, None, None)
	else: # Iterate, check for overlaps - merge or add as new.

		isOverlapFound = False
		for idx, (startSaved, endSaved, ISidSaved, orientSaved) in enumerate(dict_reducedSet[refId][SITES]):
			if ISidSaved == ISid and orientSaved == orient and isAnyOverlap(startSaved, endSaved, start_min, end_max, 0) :

				# print ("Merge: - i.e. update startSaved & endSaved")

				newStart = startSaved if (startSaved <= start_min) else start_min
				newEnd = endSaved if endSaved >= end_max else end_max


				dict_reducedSet[refId][SITES][idx] = (newStart, newEnd, ISidSaved, orientSaved)
				addSupportingInfo(dict_reducedSet, refId, newStart, newEnd, score, side_s, side_l, startSaved, endSaved, ISid, orient, ISidSaved, orientSaved)


				isOverlapFound = True
				break

		if isOverlapFound == False:
			dict_reducedSet[refId][SITES].append((start_min, end_max, ISid, orient))
			addSupportingInfo(dict_reducedSet, refId, start_min, end_max, score, side_s, side_l, None, None, ISid, orient, None, None)


def addSupportingInfo_reduced(dict_info_new, dict_info_prev, newStart, newEnd, prevStart, prevEnd, oldStart, oldEnd, newISid, newOrient, prevISid, prevOrient, oldISid, oldOrient):

	if (newStart, newEnd, newISid, newOrient) not in dict_info_new:
		dict_info_new[(newStart, newEnd, newISid, newOrient)] = {}
		dict_info_new[(newStart, newEnd, newISid, newOrient)][NUM_INF] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][NUM_INF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_L]
        
	else: 
		dict_info_new[(newStart, newEnd, newISid, newOrient)][NUM_INF] = dict_info_new[(newStart, newEnd, newISid, newOrient)][NUM_INF] +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][NUM_INF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE]  +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] +  dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] + dict_info_prev[(prevStart, prevEnd, prevISid, prevOrient)][SIDE_L]

	if oldStart != None and oldEnd != None and oldISid != None and oldOrient != None and not (newStart == oldStart and newEnd == oldEnd and newISid == oldISid and newOrient == oldOrient): 
		# print ("Now also add, this info, then del key.")

		dict_info_new[(newStart, newEnd, newISid, newOrient)][NUM_INF] = dict_info_new[(newStart, newEnd, newISid, newOrient)][NUM_INF] +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][NUM_INF]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE] = dict_info_new[(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE]  +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][COMBINED_SCORE]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_S] +  dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][SIDE_S]
		dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_info_new[(newStart, newEnd, newISid, newOrient)][SIDE_L] + dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)][SIDE_L]

		del dict_info_new[(oldStart, oldEnd, oldISid, oldOrient)]


def addSupportingInfo(dict_, contigId, newStart, newEnd, score, side_s, side_l, oldStart, oldEnd, newISid, newOrient, oldISid, oldOrient):

	if oldStart == None and oldEnd == None:
		dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)] = {NUM_INF: 0, COMBINED_SCORE: 0, SIDE_S: 0, SIDE_L: 0}

	else:
		if (oldStart, oldEnd, oldISid, oldOrient) in dict_[contigId][INFO]:

			if (oldStart == newStart and oldEnd == newEnd and oldISid == newISid and oldOrient == newOrient) :
				# print ('Simply update: i.e. below')

				pass
			else:
				# print ('Delete the keys')
				dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)] = {NUM_INF: dict_[contigId][INFO][(oldStart, oldEnd, oldISid, oldOrient)][NUM_INF], COMBINED_SCORE: dict_[contigId][INFO][(oldStart, oldEnd, oldISid, oldOrient)][COMBINED_SCORE], SIDE_S: dict_[contigId][INFO][(oldStart, oldEnd, oldISid, oldOrient)][SIDE_S], SIDE_L: dict_[contigId][INFO][(oldStart, oldEnd, oldISid, oldOrient)][SIDE_L]}

				del dict_[contigId][INFO][(oldStart, oldEnd, oldISid, oldOrient)]


		# elif (oldStart == newStart and oldEnd == newEnd)


		else:

			print (dict_[contigId][INFO].keys())
			print ('hmmm key not found... ' + str(oldStart) + ':' + str(oldEnd) + " New: " + str(newStart) + ":" + str(newEnd))


	if score == 'Inf':
		dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][NUM_INF] = dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][NUM_INF] + 1
	else:
		# print ('Score: ' + score)
		dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE] = dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][COMBINED_SCORE] + int(score)

	if side_s != None:
		dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][SIDE_S] = dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][SIDE_S] + int(side_s)

	if side_l != None:
		dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][SIDE_L] = dict_[contigId][INFO][(newStart, newEnd, newISid, newOrient)][SIDE_L] + int(side_l)



def reduceSet_(dict_allISsites, th_forMergingOverlaps):
	dict_reducedSet = dict() # dict_[refId|contigId] => (start, end, IS, orient)

	list_allISsites = sorted(dict_allISsites.keys(), key=lambda x: x[1])

	for (refId, startPos, endPos, ISid, orient, isCalcFromContig) in list_allISsites:


		print("refId: " + refId + " start: " + str(startPos) + " end: " + str(endPos))



	"""
		# isAnyOverlap()
		if len(dict_reducedSet.keys()) == 0:
			dict_reducedSet[(refId, startPos, endPos, ISid, orient)] = dict()
			dict_reducedSet[(refId, startPos, endPos, ISid, orient)]['count'] = 1
			dict_reducedSet[(refId, startPos, endPos, ISid, orient)]['isCalcFromContig'] = [isCalcFromContig]

		else: # loop thru and check
			list_reducedSites = sorted(dict_reducedSet.keys(), key=lambda x: x[0])
			isFound = False
			for (rs_start, rs_end, rs_ISid, rs_orient) in list_reducedSites:

				if rs_ISid == ISid and rs_orient == orient and isAnyOverlap(startPos, endPos, rs_start, rs_end, th_forMergingOverlaps) == True :

					isFound = True


					dict_reducedSet[(rs_start, rs_end, rs_ISid, rs_orient)]['count'] = dict_reducedSet[(rs_start, rs_end, rs_ISid, rs_orient)]['count'] + 1
					dict_reducedSet[(rs_start, rs_end, rs_ISid, rs_orient)]['isCalcFromContig'].append(isCalcFromContig)

					# print(dict_allISsites[(startPos, endPos, ISid, orient, isCalcFromContig)])


			if isFound == False:
				dict_reducedSet[(startPos, endPos, ISid, orient)] = dict()
				dict_reducedSet[(startPos, endPos, ISid, orient)]['count'] = 1
				dict_reducedSet[(startPos, endPos, ISid, orient)]['isCalcFromContig'] = [isCalcFromContig]

	for key in dict_reducedSet:
		print (str(key) + '\t' + str(dict_reducedSet[key]['count']))
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


############################################ AUX - LOAD DATA
def loadTheData(fn_ISinGff):

	COL_ID = 0
	COL_START = 3
	COL_END = 4
	COL_SCORE = 5
	COL_ORIENT = 6
	COL_INFO = 8

	COL_INFO_NAME = 0


	dict_ISinsertions = dict() # dict_[(refId|contigId, startPos, endPos, IStype, orient)] => [(modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, score)]
	dict_seqRegInfo = dict() # dict_[refId|contigId] => (start, end)
	dict_colors = dict() 

	with open(fn_ISinGff, 'r') as fh:
		for line in fh:

			if re.match(r'^\#\#gff-version', line):
				continue

			elif re.match(r'^\#\#sequence-region', line):
				extractAndAddSeqRegion(dict_seqRegInfo, line)
				# print (line)
				continue

			line = line.strip()

			arr = line.split('\t')

			# print (line)


			startPos = int(arr[COL_START])
			endPos = int(arr[COL_END])

			score = arr[COL_SCORE]
			orient = arr[COL_ORIENT]



			arr_info = arr[COL_INFO].split(';')
			# IS_type = re.sub('^Name=', '', arr_info[COL_INFO_NAME])
			# isCalcFromContig = re.sub('^isCalcFromContig=', '', arr_info[len(arr_info)-1])

			(modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, isCalcFromContig, IStype, color, orient_neg, orient_pos) = loadInfoColumn(arr_info)

			# print (str(startPos) + ":" + str(endPos) + "\t" + score + "\t" + orient + '\t' + IS_type + '\t' + isCalcFromContig + '\t' + str((modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen)))

			if color != None: 
				if IStype not in dict_colors: 
					dict_colors[IStype] = color 

			if orient != '.': 
				
				if (arr[COL_ID], startPos, endPos, IStype, orient, isCalcFromContig) not in dict_ISinsertions:
					dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient, isCalcFromContig)] = []

				dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient, isCalcFromContig)].append((modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, score))

			else: # orient = None ## Need to add this twice - orient_pos and orient_neg! 
				if (arr[COL_ID], startPos, endPos, IStype, orient_neg, isCalcFromContig) not in dict_ISinsertions:
					dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient_neg, isCalcFromContig)] = []

				dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient_neg, isCalcFromContig)].append((modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, score))

				if (arr[COL_ID], startPos, endPos, IStype, orient_pos, isCalcFromContig) not in dict_ISinsertions:
					dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient_pos, isCalcFromContig)] = []

				dict_ISinsertions[(arr[COL_ID], startPos, endPos, IStype, orient_pos, isCalcFromContig)].append((modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, score))


			# if (startPos, endPos) in dict_allData:

			#    dict_allData[int(arr)]

	return (dict_ISinsertions, dict_seqRegInfo, dict_colors)

def extractAndAddSeqRegion(dict_seqRegInfo, line):
	arr = re.split(r'[\s\t]+', line)

	COL_SEQ_ID = 1
	COL_SEQ_START = 2
	COL_SEQ_END = 3

	if arr[COL_SEQ_ID] not in dict_seqRegInfo:
		dict_seqRegInfo[arr[COL_SEQ_ID]] = (int(arr[COL_SEQ_START]), int(arr[COL_SEQ_END]))


def loadInfoColumn(arr_info):
	modePos = None
	meanPos = None
	medianPos = None
	maxPos = None
	minPos = None
	side_s = None
	side_l = None
	meanAlignLen = None
	medianAlignLen = None
	isCalcFromContig = None
	IStype = None
	color = None
	orient_pos = None 
	orient_neg = None 

	for anInfo in arr_info:
		if re.match('^modeCount=', anInfo):
			modePos = re.sub('^modeCount=', '', anInfo)
			modePos = int(modePos)
		elif re.match('^meanPos=', anInfo):
			meanPos = re.sub('^meanPos=', '', anInfo)
			meanPos = float(meanPos)
		elif re.match('^medianPos=', anInfo):
			medianPos = re.sub('^medianPos=', '', anInfo)
			medianPos = float(medianPos)
		elif re.match('^maxPos=', anInfo):
			maxPos = re.sub('^maxPos=', '', anInfo)
			maxPos = int(maxPos)
		elif re.match('^minPos=', anInfo):
			minPos = re.sub('^minPos=', '', anInfo)
			minPos = int(minPos)
		elif re.match('^side_s=', anInfo):
			side_s = re.sub('^side_s=', '', anInfo)
			side_s = int(side_s)
		elif re.match('^side_l=', anInfo):
			side_l = re.sub('^side_l=', '', anInfo)
			side_l = int(side_l)
		elif re.match('^meanAlignLen=', anInfo):
			meanAlignLen = re.sub('^meanAlignLen=', '', anInfo)
			meanAlignLen = float(meanAlignLen)
		elif re.match('^medianAlignLen=', anInfo):
			medianAlignLen = re.sub('^medianAlignLen=', '', anInfo)
			medianAlignLen = float(medianAlignLen)
		elif re.match('^isCalcFromContig=', anInfo):
			isCalcFromContig = re.sub('^isCalcFromContig=', '', anInfo)
			# isCalcFromContig = float(medianAlignLen)
		elif re.match('^Name=', anInfo):
			IStype = re.sub('^Name=', '', anInfo)
		elif re.match('^color=', anInfo): 
			color = re.sub('^color=', '', anInfo)
		elif re.match('^orient\_\-', anInfo): 
			orient_neg = re.sub('^orient\_\-=', '', anInfo)
			orient_neg = '-' # int(orient_neg) 
		elif re.match('^orient\_\+', anInfo): 
			orient_pos = re.sub('^orient\_\+=', '', anInfo)
			orient_pos = '+' # int(orient_pos) 

		else: # Handling any missing cases!
			sys.stdout.write("please add this column too! " + anInfo)
			sys.exit()

	return (modePos, meanPos, medianPos, maxPos, minPos, side_s, side_l, meanAlignLen, medianAlignLen, isCalcFromContig, IStype, color, orient_neg, orient_pos)


############################################ MAIN
def main():
	parser = argparse.ArgumentParser(description='Get merged count of final IS insertions')

	parser.add_argument('--fn_ISinGff', nargs=1, required=True, help='ISinContigs.gff file, or IStoRef.gff file.')
	parser.add_argument('--th_forMergingOverlaps', nargs=1, required=False, help='Distance between two IS, before calling them merges, Default = 20.', default=[20])
	parser.add_argument('--ignoreOrient', nargs=1, required=False, help='Values as True or False - to igore or-else check orientation, repsectively. Ignored orientation is appended as orient1:orient2. Default = False', default=["False"])
	parser.add_argument('--ignoreIStype', nargs=1, required=False, help='Values as True or False - to igore or-else check IStype, repsectively. Ignored IStypes are appended as IStype1:IStype2:..:IStypeN. This automatically sets isIgnoreOrient to "True". Default = False.', default=["False"])
	parser.add_argument('--th_toCountAsEdge', nargs=1, required=False, default=[100], help='This distance from the (start+th) or (end-th) is checked to determine if insertion is counted as being an edge insertion or an insertion in the middle. Default = 100 (bps).')
	parser.add_argument('--fnOut_estimates', nargs=1, required=True, help='File name to add the estimated number of IS to.')
	parser.add_argument('--fnOut_estimates_singleRow', nargs=1, required=True, help='File name to add the estimated number of IS to.')
	parser.add_argument('--fnOut_gff3_merged', nargs=1, required=True, help='File name to add the merged positions to.')
	parser.add_argument('--separator', required=False, default=[':'], help='Separate unresolvable IS and orientations. Default is ":". Choose a separator that is not an alphabet in the IS identifier.')

	parser.add_argument('--isolateId', nargs=1, required=True, help='Identifier of currently isolate.')



	args = parser.parse_args()
	args.th_forMergingOverlaps[0] = int(args.th_forMergingOverlaps[0])

	getTheFinalCountOfIS(args.fn_ISinGff[0], int(args.th_forMergingOverlaps[0]), args.ignoreOrient[0], int(args.th_toCountAsEdge[0]), args.ignoreIStype[0], args.fnOut_estimates[0], args.separator[0], args.fnOut_gff3_merged[0], args.isolateId[0], args.fnOut_estimates_singleRow[0])

if __name__=='__main__':
	main()
