#!/usr/bin/python3

import argparse
import re
import os
import glob
import sys

PANISA = 'Panisa/'
TEFINDER = 'Tefinder/'
ISMAPPER = 'Ismapper/'
MGEFINDER_DENOVO = 'Mgefinder_denovo/'
MGEFINDER_DATABASE = 'Mgefinder_database/'
WIIS = 'Wiis_v2/'


########################################## TOP_LVL
def createSummary(IStoRef_blastRes, IStoRef_gff, th_toMergePosFound, th_finalAlignOverlap):

	# 2. For each compl genome in map file,


	dict_ISinCompl = loadISinRef(IStoRef_blastRes)
	# print (arr[illuminaCol])

	dict_wiisWrtCompl = loadWiis(IStoRef_gff)

	(dict_wiisPosMerges, allMergIds) = mergeTheOverlaps(dict_wiisWrtCompl, th_toMergePosFound)
	isInComplGenome(dict_wiisPosMerges, dict_ISinCompl, dict_wiisWrtCompl, allMergIds, th_toMergePosFound, th_finalAlignOverlap) 


	#for (start_wiis, end_wiis, isMerged) in sorted(dict_wiisWrtCompl, key=lambda tup: tup[0]): 
	#	print (str(start_wiis) + ":" + str(end_wiis) + ':' + str(isMerged) + '\t' + str(dict_wiisWrtCompl[(start_wiis, end_wiis, isMerged)])) 



			



########################################## AUX - WIIS
def isInComplGenome(dict_wiisPosMerges, dict_ISinCompl, dict_wiisWrtCompl, allMergeIds, th_toMergePosFound, th_finalAlignOverlap): 

	list_mergedIds_found = list() 
	# dict_ISinComplFound = dict() # {(refStart, refEnd)} => mergeId
	# list_ISinComplFound = list() 

	print("ISinCompl" + '\t' + 'ISfoundByWiis\t...')
	for (refStart, refEnd) in sorted(dict_ISinCompl, key=lambda tup: int(tup[0])): 
		sys.stdout.write(str(refStart) + ":" + str(refEnd))
		for (ISid, orient) in dict_ISinCompl[(refStart, refEnd)]: 
			sys.stdout.write('|' + ISid + "_" + orient)
		

		isFound = False
		for (wiisStart, wiisEnd) in dict_wiisPosMerges: 
			# if merged_id not in list_mergedIds_foundInRef:

			if isOverlap(refStart, refEnd, wiisStart, wiisEnd, th_finalAlignOverlap): 

				isFound = True 

				if dict_wiisPosMerges[(wiisStart, wiisEnd)] not in list_mergedIds_found: 
					list_mergedIds_found.append(dict_wiisPosMerges[(wiisStart, wiisEnd)])

				list_ISids = theISids(dict_wiisWrtCompl, wiisStart, wiisEnd)
				sys.stdout.write('\t' + str(wiisStart) + ":" + str(wiisEnd))
				for (ISid, orient) in list_ISids: 
					sys.stdout.write('|' + ISid + "_" + orient)
					# dict_ISinComplFound[(refStart, refEnd)].append(merged_id) 



		sys.stdout.write('\n')



	print("# Total mergeIds are (IS within " + str(th_toMergePosFound) + " bps):" + str(len(allMergeIds)))
	print("# Total mergeIds found in the complGenome: " + str(len(list_mergedIds_found)))

	list_mergedIds_notFound = set(allMergeIds).difference(set(list_mergedIds_found))
	print("# Total mergeIds not found in complGenome: " + str(len(list_mergedIds_notFound)))

	print("# The mergeIds not found are: ")
	for (start, end) in dict_wiisPosMerges: 
		if (dict_wiisPosMerges[(start, end)] in list_mergedIds_notFound): 
			sys.stdout.write (str(start) + ":" + str(end)) 
			for (ISid, orient) in dict_wiisWrtCompl[(start, end)]: 
				sys.stdout.write("|" + ISid + '_' + orient)

			sys.stdout.write('\n')



def theISids(dict_wiisWrtCompl, start, end):
	return (dict_wiisWrtCompl[(start, end)])



def mergeTheOverlaps(dict_wiisWrtCompl, th_toMergePosFound): 
	allMergeIds = [] 

	theKeys = sorted(dict_wiisWrtCompl, key=lambda tup: int(tup[0])) # list{(start, end, isMerged)}
	dict_wiisPosMerges = dict() # dict_[(start, end)] => mergeId

	
	for idx in range(0, len(theKeys)): 

		currIdx = idx 
		if theKeys[idx] not in dict_wiisPosMerges: 
			dict_wiisPosMerges[theKeys[idx]] = currIdx
			allMergeIds.append(idx)

		else:
			currIdx = dict_wiisPosMerges[theKeys[idx]]

		for idx_secondloop in range(idx+1, len(theKeys)): 


			if theKeys[idx_secondloop] not in dict_wiisPosMerges: 

				if isOverlap(theKeys[idx][0], theKeys[idx][1], theKeys[idx_secondloop][0], theKeys[idx_secondloop][1], th_toMergePosFound): 
					dict_wiisPosMerges[theKeys[idx_secondloop]] = currIdx
					

	# s = set(dict_wiisPosMerges.values())
	#for val in dict_wiisPosMerges.values():
	#	s.add(val)

	# print('Uniq mergeIds are ' + str(len(allMergeIds)))
	# print('Count in dict_values ' + str(len(s)))



	return (dict_wiisPosMerges, allMergeIds)



def isOverlap(s1, e1, s2, e2, th): 
	if (s1 >= (s2 - th) and s1 <= (e2 + th)) or (e1 >= (s2 - th) and e1 <= (e2 + th)) or (s2 >= (s1 - th) and s2 <= (e1 + th)) or (e2 >= (s1 - th) and e2 <= (e1 + th)): 
		return True 

	return False 

########################################## AUX - WIIS
def loadWiis(IStoRef_gff):

	dict_wiis = {} # dict_{(refStart,refEnd)} = [(ISid, orient, isFoundInRef)]

	fn = IStoRef_gff

	if not os.path.exists(fn):
		return dict_wiis


	Col_start = 3
	Col_end = 4
	Col_orient = 6
	Col_info = 8
	Col_score = 5


	with open(fn, 'r') as fh:
		for line in fh:

			if re.match("^\#", line):
				continue

			
			line = line.strip()
			arr = re.split("\t", line)

			if 'isCalcFromContig=True' in arr[Col_info]:
				
				# print (line)
				ISid = ''
				arr_info = arr[Col_info].split(';')
				for val in arr_info:
					# print (arr_info)
					if re.match('^Name', val):
						ISid = val

				ISid = re.sub('^Name\=', '', ISid)

				# print (arr[Col_start] + ':' + arr[Col_end] + '\t' + arr[Col_orient] + '\t' + ISid)

				startPos = int(arr[Col_start])
				endPos = int(arr[Col_end]) 

				if (startPos, endPos) not in dict_wiis:
					dict_wiis[(startPos, endPos)] = []

				if (ISid, arr[Col_orient]) not in dict_wiis[(startPos, endPos)]: 
					dict_wiis[(startPos, endPos)].append((ISid, arr[Col_orient]))

	# for key in dict_wiis:
	# 	print (str(key) + '\t' + str(dict_wiis[key]))
	return dict_wiis


########################################## AUX - MGEFINDER


########################################## AUX - TEFINDER

########################################## AUX - PANISA
def findAndPrint_panisa(dict_ISinRef, dict_panisa, dict_toolSeen, fh_out):

	total_foundIS = 0
	total_foundSameIS = 0
	total_foundSameOrient = 0

	ISinRefPositions = list(dict_ISinRef.keys())
	ISinRefPositions_sorted = sorted(ISinRefPositions, key=lambda x: x[0])

	# print(ISinRefPositions_sorted)
	for (refStart, refEnd) in ISinRefPositions_sorted:
		for (ISid, orient) in dict_ISinRef[(refStart, refEnd)]:
			(isISfound, isSameIS, isSameOrient) = isInTool(dict_panisa, refStart, refEnd, ISid, orient, dict_toolSeen)
			# print (str(isFound) + '\t' + str(isSameIS) + '\t' + str(isSameOrient))

			if isISfound == True:
				total_foundIS = total_foundIS + 1

				if isSameIS == True:
					total_foundSameIS = total_foundSameIS + 1

					if isSameOrient == True:
						total_foundSameOrient = total_foundSameOrient + 1




			printToOutfile_tsv(fh_out, refStart, refEnd, ISid, orient, isISfound, isSameIS, isSameOrient)


	fh_out.write('\t' + str(total_foundIS) + '\t' + str(total_foundSameIS) + '\t' + str(total_foundSameOrient))


	totalNotSeen = calcNotSeen(dict_panisa, dict_toolSeen)

	fh_out.write('\t' + str(totalNotSeen))
	#for (panisaStart, panisaEnd) in dict_panisa:
	#	isFound = isInRef(dict_ISinRef, panisaStart, panisaEnd)
	#	print (isFound)

	#print ('Printing the dict_toolsSeen ' + str(len(dict_toolSeen.keys())))
	#for key in dict_toolSeen:
	#	print (str(key) + '\t' + str(dict_toolSeen[key]))


def loadPanISa(SRRid, dir_panISa, dir_runs):

	dict_panisa = {} # dict_{(refStart, refEnd)} = [(ISid, orient, isFoundInRef)]

	fn_panISaRes = dir_runs + SRRid + '/' + dir_panISa + 'panISa_withISFinder_output'

	Col_refStart = 2
	Col_refEnd = 3
	Col_ISid = 5

	if (os.path.isfile(fn_panISaRes)):
		with open (fn_panISaRes, 'r') as fh:

			isHeader = True
			for line in fh:


				if isHeader == True:
					isHeader = False
					continue

				line = line.strip()

				arr = line.split("\t")

				orient = calcOrient_panISa(int(arr[Col_refStart]), int(arr[Col_refEnd]))

				# print (line)
				# print (arr[Col_refStart] + ":" + arr[Col_refEnd] + '\t' + arr[Col_ISid] + '\t' + orient)

				refStart = int(arr[Col_refStart])
				refEnd = int(arr[Col_refEnd])

				if (refStart, refEnd) not in dict_panisa:
					dict_panisa[(refStart, refEnd)] = []

				dict_panisa[(refStart, refEnd)].append((arr[Col_ISid], orient))

	return dict_panisa

def calcNotSeen(dict_tool, dict_toolSeen):

	totalNotSeen = 0

	for (start, end) in dict_tool:
		if (start, end) not in dict_toolSeen:
			totalNotSeen = totalNotSeen + 1

	return totalNotSeen

########################################## AUX - common for all tools

def printToOutfile_tsv(fh_out, refStart, refEnd, ISid, orient, isISfound, isSameIS, isSameOrient):

	if isISfound:
		fh_out.write("\t1")

		if isSameIS:
			fh_out.write("_1")

			if isSameOrient:
				fh_out.write("_1")

	else:
		fh_out.write("\t0")


def isInTool(dict_tool, refStart, refEnd, refISid, refOrient, dict_toolSeen):
	isISfound = False; isSameIS = False; isSameOrient = False;

	for (toolStart, toolEnd) in dict_tool:
		if (toolStart >= refStart and toolStart <= refEnd) or (toolEnd >= refStart and toolEnd <= refEnd) or (refStart >= toolStart and refStart <= toolEnd) or (refEnd >= toolStart and refEnd <= toolEnd):
			isISfound = True

			# Adding to tool
			if (toolStart, toolEnd) not in dict_toolSeen:
				dict_toolSeen[(toolStart, toolEnd)] = []


			for (toolISid, toolOrient) in dict_tool[(toolStart, toolEnd)]:
				if re.match(toolISid, refISid, flags=re.I) or re.match(refISid, toolISid, flags=re.I):
					isSameIS = True

					# Adding to tool
					if (toolISid, toolOrient) not in dict_toolSeen[(toolStart, toolEnd)]:
						dict_toolSeen[(toolStart, toolEnd)].append((toolISid, toolOrient))

				if toolOrient == refOrient:
					isSameOrient = True


					# Adding to tool
					if (toolISid, toolOrient) not in dict_toolSeen[(toolStart, toolEnd)]:
						dict_toolSeen[(toolStart, toolEnd)].append((toolISid, toolOrient))
			break

	return (isISfound, isSameIS, isSameOrient)


########################################## AUX
def printHeader(fh_out, dict_ISinRef):
	ISinRefPositions = list(dict_ISinRef.keys())
	ISinRefPositions_sorted = sorted(ISinRefPositions, key=lambda x: x[0])

	fh_out.write('Reference')

	for (refStart, refEnd) in ISinRefPositions_sorted:
		for (ISid, orient) in dict_ISinRef[(refStart, refEnd)]:
			fh_out.write('\t' + str(refStart) + ":" + str(refEnd) + "_" + ISid + "_" + orient)

	fh_out.write('\t' + 'TotalISfound' + '\t' + 'OfWhichSameIS' + '\t' + 'OfWhichSameOrient' + '\t' + 'ToolISposNotInRef' '\n')


def calcOrient_panISa(refStart, refEnd):

	if (refStart > refEnd):
		return '-'

	return '+'

def loadISinRef(fn_IStoRef):
	dict_ISinRef = {} # [(start, end)] => [(ISid, orient), (ISid, orient)..[ ]

	Col_startInRef = 8
	Col_endInRef = 9
	Col_ISid = 0

	Col_startOfIS = 6
	Col_endOfIS = 7

	with open(fn_IStoRef, 'r') as fh:
		for line in fh:
			if not re.match('^\#', line):

				line = line.strip()

				arr = line.split('\t')

				refStart = int(arr[Col_startInRef])
				refEnd = int(arr[Col_endInRef])


				orient = calcOrient(refStart, refEnd, int(arr[Col_startOfIS]), int(arr[Col_endOfIS]))


				if refStart > refEnd:
					tmp = refStart
					refStart = refEnd
					refEnd = tmp


				# Add to dict:
				if (refStart, refEnd) not in dict_ISinRef:
					dict_ISinRef[(refStart, refEnd)] = []

				dict_ISinRef[(refStart, refEnd)].append((arr[Col_ISid], orient))


				# print (arr[Col_startInRef] + ':' + arr[Col_endInRef] + '\t' + arr[Col_ISid] + '\t' + arr[Col_startOfIS] + ":" + arr[Col_endOfIS] + '\t' + orient + "\t" + str(refStart) + ":" + str(refEnd))

	return dict_ISinRef

def calcOrient(refStart, refEnd, ISstart, ISend):

	if (refStart >= refEnd and ISstart >= ISend) or (refStart <= refEnd and ISstart <= ISend):
		return '+'

	return '-'

def checkAndAddSlash(dir_out):

	if not re.search(r'/$', dir_out):
		dir_out = dir_out + "/"

	return dir_out

########################################## TOP_LVL

def main():
	parser = argparse.ArgumentParser(description='Count IS found by WiIS vs the complete genome.')

	# parser.add_argument('--mapFile', nargs=1, required=True, help='Map file containing the SRR to the complete genomes')
	# parser.add_argument('--map_illuminaCol', nargs =1, required=True, type=int, help='Column number of Illumina reads name')
	# parser.add_argument('--map_complGenomeCol', nargs =1, required=True, type=int, help='Column number of complete genome name')
	# parser.add_argument('--dir_runs', nargs=1, required=True, help='Output folder for run outputs.')


	parser.add_argument('--IStoRef_blastRes', nargs=1, required=True, help='E.g. CompleteGenome/IStoComple.blastRes')

	parser.add_argument('--IStoRef_gff', nargs=1, required=True, help='E.g. /pathToIStoRef/IStoRef.gff')

	# parser.add_argument('--outputFile', nargs=1, required=True, help='E.g. wiis_compl.txt')

	parser.add_argument('--th_toMergePosFound', nargs=1, required=False, default=[20], help='To merge IS-positions found within \'th\'. Default=20.')

	parser.add_argument('--th_finalAlignOverlap', nargs=1, required=False, default=[0], help='To merge IS-positions found within \'th\'. Default=0.')


	# parser.add_argument('--mgefinder_denovo', nargs=1, required=False, help='Filename path to 01....tsv')
	
	# parser.add_argument('--mgefinder_database', nargs=1, required=False, help='Filename path to 01....tsv')

	args = parser.parse_args()

	# args.programName[0] = checkAndAddSlash(args.programName[0])
	# args.dir_runs[0] = checkAndAddSlash(args.dir_runs[0])


	if (not os.path.exists(args.IStoRef_blastRes[0])): 
		sys.exit("Error: " + args.IStoRef_blastRes[0] + ' not found.')

	if (not os.path.exists(args.IStoRef_gff[0])): 
		sys.exit("Error: " + args.IStoRef_gff[0] + ' not found.')	

	createSummary(args.IStoRef_blastRes[0], args.IStoRef_gff[0], int(args.th_toMergePosFound[0]), int(args.th_finalAlignOverlap[0]))


if __name__ == '__main__':
	main()
