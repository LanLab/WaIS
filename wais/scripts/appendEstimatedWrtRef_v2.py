
import argparse 
import sys
import re

################################################ TOP_LVL 
def determineAndCalcEstimated(fn_IStoRef, fn_ISinContigs, fn_estimates, fn_estimates_singleRow): 

	(list_inRef_contigISids, dict_counts_isNew, list_inRef_ISids) = loadCounts(fn_IStoRef) # Load counts of isNew=False, isNew=True (here, either simple counts, or ISid & orient), and contigIds. 
	dict_isNewIsFalse = breakFalseCountsByIS(dict_counts_isNew)

	(list_allContigIds, dict_contigInfo) = loadContigISids(fn_ISinContigs)

	# getContigIdsNotMapped()

	for refId in dict_counts_isNew: 
		for isNew in dict_counts_isNew[refId]: 
			if isNew == 'False': 
				print('Total same: ' + str(len(dict_counts_isNew[refId][isNew])))
			# dict_counts_isNew[refIS]

		for ISid in dict_isNewIsFalse[refId]: 
			print (ISid + ' ' + str(dict_isNewIsFalse[refId][ISid]))

	list_missed = getMissedContigs(list_allContigIds, list_inRef_contigISids)

	(dict_missedStats, list_posInContig_encountered) = calcMissedCounts(list_missed, dict_contigInfo)

	print (list_missed) 

	for ISid in dict_missedStats:
		for posInContig in dict_missedStats[ISid]:  
			print(ISid + ' ' + posInContig + ' ' + str(dict_missedStats[ISid][posInContig])) 

	# Printing final outputs
	(estSingleRow_header, estSingleRow_vals) = appendToEstimates(fn_estimates, list_inRef_ISids, dict_isNewIsFalse, dict_counts_isNew, list_posInContig_encountered, dict_missedStats)
	
	appendToEstimates_singleRow(fn_estimates_singleRow, estSingleRow_header, estSingleRow_vals)

################################################ OUTPUT
def appendToEstimates_singleRow(fn_estimates_singleRow, estSingleRow_header, estSingleRow_vals): 
	lines = None 

	with open(fn_estimates_singleRow, 'r') as fh:
		lines = fh.read() 

	with open (fn_estimates_singleRow, 'w+') as fh: 
		for idx, line in enumerate(lines.split('\n')): 
			if idx == 0: 
				line = line + '\t'.join(estSingleRow_header) 
			if idx == 1: 
				line = line + '\t'.join(estSingleRow_vals) 

			fh.write(line + '\n')	


def appendToEstimates(fn_estimates, list_allISids, dict_isNewIsFalse, dict_counts_isNew, list_posInContig_encountered, dict_missedStats): 
	estSingleRow_header = []  
	estSingleRow_vals = []


	# print (list_inRef_ISids)
	with open(fn_estimates, 'a+') as fh: 
		printHandle = fh # sys.stdout

	
		printHandle.write('#################### IS mapped from contigs to reference\n')
		estSingleRow_header.append('#Reference')
		estSingleRow_vals.append('') 

		printHandle.write('RefId')
		for ISid in sorted(list_allISids): 
			printHandle.write('\t' + ISid)
		printHandle.write('\tTotal\n')


		for refId in dict_isNewIsFalse: 
			
			# 1. isNew = False 
			isNew = 'False'
			total_row = 0
	
			printHandle.write(refId + ':' + "Ref_IS_found") 

			estSingleRow_header.append(refId)
			estSingleRow_vals.append('')
				
			for ISid in sorted(list_allISids): 

				printHandle.write('\t')

				estSingleRow_header.append(ISid + ':Ref_IS_found' )
				singleRow_valToAppend = '' 

				if isNew == 'False': 
					# use counts stored in dict_isNewIsFalse
					if ISid in dict_isNewIsFalse[refId]: 
						printHandle.write(str(dict_isNewIsFalse[refId][ISid]))
						total_row = total_row + dict_isNewIsFalse[refId][ISid]

						singleRow_valToAppend = dict_isNewIsFalse[refId][ISid]

				estSingleRow_vals.append(str(singleRow_valToAppend))
		

			"""		
			if ISid in dict_counts_isNew[refId][isNew]:
					printHandle.write(ISid + ':' + str(dict_counts_isNew[refId][isNew][ISid]))
			"""

			estSingleRow_header.append('Total:Ref_IS_found')
			estSingleRow_vals.append(str(total_row))
			
			printHandle.write('\t' + str(total_row))
			printHandle.write('\n')


			# 2. isNew = True 
			isNew = "True"
			total_row = 0 
			
			# for posInContig in list_posInContig_encountered: 
			printHandle.write(refId + ':' + "New_IS") #  + ':Mapped_from_' + posInContig) 


			for ISid in sorted(list_allISids): 
				printHandle.write('\t')

				estSingleRow_header.append(ISid + ':' + "New_IS")
				singleRow_valToAppend = '' 

				if refId in dict_counts_isNew and isNew in dict_counts_isNew[refId] and ISid in dict_counts_isNew[refId][isNew]:
					printHandle.write(str(dict_counts_isNew[refId][isNew][ISid]))
					total_row = total_row + dict_counts_isNew[refId][isNew][ISid]

					# estSingleRow_header.append(ISid + ':' + "New_IS")
					# estSingleRow_vals.append(dict_counts_isNew[refId][isNew][ISid])
					singleRow_valToAppend = dict_counts_isNew[refId][isNew][ISid]

				estSingleRow_vals.append(str(singleRow_valToAppend))
				

			printHandle.write('\t' + str(total_row))
			printHandle.write('\n')

			estSingleRow_header.append('Total:New_IS')
			estSingleRow_vals.append(str(total_row))

			#estSingleRow_header.append('Total' + ':' + "New_IS")
			#estSingleRow_vals.append(total_row)

			# 3. Total of isNew=True & isNew=False (Per reference :-)
			overallTotal = 0
			
			printHandle.write(refId + ':' + "Total") 
			for ISid in sorted(list_allISids): 
				printHandle.write('\t')

				estSingleRow_header.append(ISid + ':Total')


				total_ISid = 0 
				if ISid in dict_isNewIsFalse[refId]: 
					total_ISid = total_ISid + dict_isNewIsFalse[refId][ISid]

				# for posInContig in list_posInContig_encountered: 
				if refId in dict_counts_isNew and isNew in dict_counts_isNew[refId] and ISid in dict_counts_isNew[refId][isNew]:
					#if posInContig == 'Edge': 
					#	total_ISid = total_ISid + (dict_counts_isNew[refId][isNew][ISid])/2
					#else: 
					total_ISid = total_ISid + (dict_counts_isNew[refId][isNew][ISid])


				printHandle.write(str(total_ISid))	
				overallTotal = overallTotal + total_ISid
				
				estSingleRow_vals.append(str(total_ISid))

			printHandle.write('\t' + str(overallTotal) + '\n')

			estSingleRow_header.append('Total')
			estSingleRow_vals.append(str(overallTotal))
			

		# 4. New ISids not found in the reference
		printHandle.write('#################### IS not mapped from contigs to reference\n')
 
		estSingleRow_header.append('#IS_not_mapped_to_ref')
		estSingleRow_vals.append('')

		for ISid in sorted(list_allISids): 
			printHandle.write('\t' + ISid)
		printHandle.write('\t' + 'Total')
		printHandle.write('\n')


		overallTotal = 0 
		for posInContig in list_posInContig_encountered: 
			
			total_row = 0
			
			printHandle.write(posInContig)

			for ISid in sorted(list_allISids): 
				printHandle.write('\t') 

				estSingleRow_header.append(ISid + ':' + posInContig)
				singleRow_valToAppend = '' 
				
				# if refId in dict_counts_isNew and isNew in dict_counts_isNew[refId] and ISid in dict_counts_isNew[refId][isNew]:
				if ISid in dict_missedStats and posInContig in dict_missedStats[ISid]: 
					printHandle.write(str(dict_missedStats[ISid][posInContig]))
					total_row = total_row + dict_missedStats[ISid][posInContig]
					overallTotal = overallTotal + dict_missedStats[ISid][posInContig]
					singleRow_valToAppend = dict_missedStats[ISid][posInContig]

				
				estSingleRow_vals.append(str(singleRow_valToAppend))


			printHandle.write('\t' + str(total_row))		

			printHandle.write('\n')
		
		# 4a. Total of Middle-&-Edge.
		printHandle.write('Total')
		overallTotal = 0 
		for ISid in sorted(list_allISids): 
			
			total_ISid = 0 
			printHandle.write('\t')
			for posInContig in list_posInContig_encountered:
				if ISid in dict_missedStats and posInContig in dict_missedStats[ISid]: 
					if posInContig == 'Edge': 
						total_ISid = total_ISid + (dict_missedStats[ISid][posInContig])/2
					else: 
						total_ISid = total_ISid + dict_missedStats[ISid][posInContig]
						

			printHandle.write(str(total_ISid))
			overallTotal = overallTotal + total_ISid

			estSingleRow_header.append(ISid + ':' + 'Total')
			estSingleRow_vals.append(str(total_ISid))

		printHandle.write('\t' + str(overallTotal))
		printHandle.write('\n')

		estSingleRow_header.append('OverallTotal')
		estSingleRow_vals.append(str(overallTotal))



	########################### 
	print ('\t'.join(estSingleRow_header))
	print ('\t'.join(estSingleRow_vals)) 

	return (estSingleRow_header, estSingleRow_vals)

################################################ AUX
def calcMissedCounts(list_missed, dict_contigInfo): 

	dict_missedStats = dict() # dict_[ISid] => [EDGE|MIDDLE] => cnt
	list_posInContig_encounttered = list() # list(edge, middle, ...)

	for contigISid in list_missed: 
		(ISid, posInContig) = dict_contigInfo[contigISid]

		if ISid not in dict_missedStats: 
			dict_missedStats[ISid] = dict() 

		if posInContig not in dict_missedStats[ISid]: 
			dict_missedStats[ISid][posInContig] = 0 

		dict_missedStats[ISid][posInContig] = dict_missedStats[ISid][posInContig] + 1 

		if posInContig not in list_posInContig_encounttered: 
			list_posInContig_encounttered.append(posInContig) 


	return (dict_missedStats, list_posInContig_encounttered)

def getMissedContigs(list_all, list_inRef):

	list_missed = list(set(list_all) - set(list_inRef))

	# print (list_inRef)
	# print (len(list_missed))
	# print (list_missed) 

	return list_missed


def loadContigISids(fn_ISinContigs):

	list_allContigIds = list() # [contigISid1, contigISid2, ...]
	dict_contigInfo = dict() # dict_{contigISid1} => (ISid, EDGE|MIDDLE)

	Col_info = 8 

	with open(fn_ISinContigs, 'r') as fh: 
		for line in fh: 
			if not re.match(r'^\#', line): 
				line = line.strip() 

				arr = line.split('\t') 

				(uniqId, IStype, positionInContig) = extractAttrInfo_contigIds(arr[Col_info])

				# Note all contigIds 
				if uniqId not in list_allContigIds: 
					list_allContigIds.append(uniqId) 
				
				addToDict_contigISinfo(dict_contigInfo, uniqId, IStype, positionInContig)


				# print(uniqIds)
	return (list_allContigIds, dict_contigInfo)


def addToDict_contigISinfo(dict_contigInfo, uniqId, IStype, positionInContig): 

	if uniqId not in dict_contigInfo: 
		dict_contigInfo[uniqId] = (IStype, positionInContig) 
		

def extractAttrInfo_contigIds(strInfo): 
	arr_info = strInfo.split(';')

	uniqId = '' 
	IStype = '' 
	positionInContig = '' 

	for val in arr_info: 
		if re.match(r'^uniqId\=', val): 
			uniqId = re.sub('^uniqId\=', '', val)
		if re.match(r'^Name\=', val): 
			IStype = re.sub('^Name\=', '', val)
		if re.match(r'^positionInContig\=', val): 
			positionInContig = re.sub(r'^positionInContig\=', '', val)

	if uniqId == '' or IStype == '' or positionInContig == '':  
		sys.exit('Error: could not find uniqId=, Name=, or positionInContig= in --ISinContigs file.')

	return (uniqId, IStype, positionInContig)

def breakFalseCountsByIS(dict_counts_isNew): 

	dict_isNewIsFalse = dict() # dict_{refId} => {ISid} = cnt

	for refId in dict_counts_isNew: 

		if refId not in dict_isNewIsFalse: 
			dict_isNewIsFalse[refId] = dict() 
	
		if 'False' in dict_counts_isNew[refId]:  
			for refIS in dict_counts_isNew[refId]['False']: 
				(start, end, ISid, orient) = refIS.split(':')

				if ISid not in dict_isNewIsFalse[refId]: 
					dict_isNewIsFalse[refId][ISid] = 0 

				dict_isNewIsFalse[refId][ISid] = dict_isNewIsFalse[refId][ISid] + 1 
				

	return dict_isNewIsFalse



def loadCounts(fn_IStoRef): 

	dict_counts_isNew = dict() # dict_[refId] => {True} => {ISid} => cnt 
	# if isNew=False: dict_[refId] => {False} => [refIS1, refIS2, ... ]

	list_contigISids = list() # [contigISid1, contigISid2, contigISid3 .... ]

	list_inRef_ISids = list() # [ISid1, ISid2, ... ISidN]

	Col_refId = 0
	Col_arrInfo = 8 

	with open(fn_IStoRef, 'r') as fh: 
		for line in fh: 

			if re.match('^\#', line): 
				continue 
			
			line = line.strip()

			arr = line.split('\t') 

			(name, refIS, directInContig, isNew) = extractAttrInfo(arr[Col_arrInfo]) 


			# Adding to dict_ and list_
			addTheDataToList(list_contigISids, directInContig)
			addTheDataToDict(dict_counts_isNew, arr[Col_refId], name, refIS, isNew)
			
			addToEncounteredIS(list_inRef_ISids, isNew, refIS, name) 
			# print (line) 

	return (list_contigISids, dict_counts_isNew, list_inRef_ISids)


def addToEncounteredIS(list_inRef_ISids, isNew, refIS, ISid): 

	if isNew == 'False': 
		arr_refIS = refIS.split(',')
		for val in arr_refIS: 
			(start, end, ISid, orient) = val.split(':')

			if ISid not in list_inRef_ISids: 
				list_inRef_ISids.append(ISid)

	if isNew == 'True': 
		if ISid not in list_inRef_ISids: 
			list_inRef_ISids.append(ISid) 



def addTheDataToList(list_toAppendTo, strContigISid): 
	arr = strContigISid.split(',')

	for val in arr: 
		if val not in list_toAppendTo: 
			list_toAppendTo.append(val)



def addTheDataToDict(dict_toAddTo, refId, ISid, refIS, isNew): 

	if refId not in dict_toAddTo: 
		dict_toAddTo[refId] = dict() 

	if isNew not in dict_toAddTo[refId]: 
		if isNew == 'False': 
			dict_toAddTo[refId][isNew] = list()

		if isNew == 'True': 
			dict_toAddTo[refId][isNew] = dict()

	
	if isNew == 'False': 
		if refIS not in dict_toAddTo[refId][isNew]: 
			arr_refIS = refIS.split(',')	

			for val in arr_refIS: 
				if val not in dict_toAddTo[refId][isNew]: 
					dict_toAddTo[refId][isNew].append(val) 


	if isNew == 'True': 
		if ISid not in dict_toAddTo[refId][isNew]: 
			dict_toAddTo[refId][isNew][ISid] = 0 

		dict_toAddTo[refId][isNew][ISid] = dict_toAddTo[refId][isNew][ISid] + 1
			 

	

def extractAttrInfo(strInfo): 

	name = '' 
	refIS = '' 
	ISids_inContigs = ''
	isNew = '' 

	arr = strInfo.split(';')

	for val in arr: 
		if re.match('^Name=', val): 
			name = re.sub('^Name=', '', val)
		if re.match('isNew=', val): 
			isNew = re.sub('^isNew=', '', val) 
		if re.match('^refIS=', val): 
			refIS = re.sub('^refIS=', '', val)
		if re.match('^ISids_inContigs=', val): 
			ISids_inContigs = re.sub('^ISids_inContigs=', '', val)			

	if name == '' or isNew == '' or ISids_inContigs == '': 
		sys.exit('Error: in file --IStoRef either Name=, isNew= or DirectInContig= is missing')

	if isNew == 'False' and refIS == '': 
		sys.exit('Error: in file --IStoRef either when isNew==False, then refIS is missing.')

	# print (ISids_inContigs)
	return (name, refIS, ISids_inContigs, isNew) 


################################################ MAIN
def main(): 
	parser = argparse.ArgumentParser(description='Extract estimates w.r.t. reference and append to assembly estimates file.')

	parser.add_argument('--IStoRef', nargs=1, required=True, help='ISmappedToRef_.gff3 (with attribute isNew=T|F) and either _merged, _ignoreOrient or _ignoreIStype')
	parser.add_argument('--ISinContigs', nargs=1, required=True, help='ISinContigs_.gff3 (with uniqId) to pick up those which may be missed in the mapping to ref.')
	
	parser.add_argument('--fn_estimates', nargs=1, required=True, help="Estimates file to append the estimated counts to.")
	parser.add_argument('--fn_estimates_singleRow', nargs=1, required=True, help="Estimates file (in single row format) to append the estimated counts to.")


	# parser.add_argument('--th_overlap', nargs=1, default=[20], type=int, help='Overlaps between ISids to consider as one. Default=20.')
	# parser.add_argument('--separator', default=[':'])
	## parser.add_argument('--isMerged', action='store_true', help='Boolean specifying if ISmappedToRef_.gff3 is calculated with merged.')
	##3 parser.add_argument('--isIgnoreOrient', action='store_true', help='Boolean specifying if ISmappedToRef_.gff3 is calculated with ignoreOrient.')
	## parser.add_argument('--isIgnoreIStype', action='store_true', help='Boolean specifying if ISmappedToRef_.gff3 is calculated with ignoreIStype.')

	args = parser.parse_args() 

	# if int(args.isMerged) + int(args.isIgnoreOrient) + int(args.isIgnoreIStype) != 1: 
	# 	sys.exit("Error: Please provide exactly one of --isMerged, --isIgnoreOrient or --isIgnoreIStype.") 
	
	determineAndCalcEstimated(args.IStoRef[0], args.ISinContigs[0], args.fn_estimates[0], args.fn_estimates_singleRow[0]) 

if __name__=='__main__': 
	main()