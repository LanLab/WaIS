#!/usr/bin/python3 

import argparse
import sys
import os 
import re
import subprocess 
import logging 
from datetime import datetime
import glob
import shutil 

logFile = 'wais_' + str(datetime.now()) + '.log'
logging.basicConfig(filename=logFile, level=logging.DEBUG)

class Thresholds:

	dict_th = {
		'getFlankingSeqs': {'--overlap': None, '--minChoppedLen': None},
		'rmFlanks_whenOneDirFullAlign_v2': {'--minChoppedLen': None, '--buffer': None}, 
		'summarizeOnlyToContig': {'--th_minPident': None, '--th_alignAndLenDiff': None},
		'calcInContig_posISOrient': {'--th_minPident': None, '--th_minPalignLen': None, '--th_minAlignLen': None, '--th_minPident_direct': None, '--th_minAlignLen_direct': None, '--kmeans_clus_start': None, '--kmeans_clus_end': None, '--th_minFlankDepth': None},
		'mergeLocalCounts': {'--th_forMergingOverlaps': None, '--th_toCountAsEdge': None}, 
		'calcInRef_posISorient_v2': {'--th_alignDiff_IStoContig': None, '--merge_th': None}, 
		'ISinRefGenome_conglomerate': {'--th_toMergePosFound': None, '--th_finalAlignOverlap': None},
		'appendEstimatedWrtRef': {'--th_overlap': None}
	}

	def __init__(self, args): 
		self.dict_th['getFlankingSeqs']['--overlap'] = args.th_getFlankingSeqs_overlap[0]
		self.dict_th['getFlankingSeqs']['--minChoppedLen'] = args.th_getFlankingSeqs_minChoppedLen[0]
		
		self.dict_th['rmFlanks_whenOneDirFullAlign_v2']['--minChoppedLen'] = args.th__rmFlanks_whenOneDirFullAlign_v2__minChoppedLen[0]
		self.dict_th['rmFlanks_whenOneDirFullAlign_v2']['--buffer'] = args.th__rmFlanks_whenOneDirFullAlign_v2__buffer[0]
		
		self.dict_th['summarizeOnlyToContig']['--th_minPident'] = args.th_summarizeOnlyToContig_minPident[0]
		self.dict_th['summarizeOnlyToContig']['--th_alignAndLenDiff'] = args.th_summarizeOnlyToContig_alignAndLenDiff[0]

		self.dict_th['calcInContig_posISOrient']['--th_minPident'] = args.th__calcInContig_posISOrient__minPident[0]
		self.dict_th['calcInContig_posISOrient']['--th_minPalignLen'] = args.th__calcInContig_posISOrient__minPalignLen[0]
		self.dict_th['calcInContig_posISOrient']['--th_minAlignLen'] = args.th__calcInContig_posISOrient__minAlignLen[0]
		self.dict_th['calcInContig_posISOrient']['--th_minPident_direct'] = args.th__calcInContig_posISOrient__minPident_direct[0]
		self.dict_th['calcInContig_posISOrient']['--th_minAlignLen_direct'] = args.th__calcInContig_posISOrient__minAlignLen_direct[0]
		self.dict_th['calcInContig_posISOrient']['--kmeans_clus_start'] = args.th__calcInContig_posISOrient__clus_start[0]
		self.dict_th['calcInContig_posISOrient']['--kmeans_clus_end'] = args.th__calcInContig_posISOrient__clus_end[0]
		self.dict_th['calcInContig_posISOrient']['--th_minFlankDepth'] = args.th__calcInContig_posISOrient__minFlankDepth[0]
		
		self.dict_th['mergeLocalCounts']['--th_forMergingOverlaps'] = args.th_mergeLocalCounts_forMergingOverlaps[0]
		self.dict_th['mergeLocalCounts']['--th_toCountAsEdge'] = args.th_mergeLocalCounts_toCountAsEdge[0]
		self.dict_th['mergeLocalCounts']['--separator'] = args.th_mergeLocalCounts_separator[0]
		
		self.dict_th['calcInRef_posISorient_v2']['--th_alignDiff_IStoContig'] = args.th__calcInRef_posISorient_v2__alignDiff_IStoContig[0]
		self.dict_th['calcInRef_posISorient_v2']['--merge_th'] = args.th__calcInRef_posISorient_v2__merge[0]
		self.dict_th['calcInRef_posISorient_v2']['--separator'] = args.th__calcInRef_posISorient_v2__separator[0]

		self.dict_th['ISinRefGenome_conglomerate']['--th_toMergePosFound'] = args.th__ISinRefGenome_conglomerate__toMergePosFound[0]
		self.dict_th['ISinRefGenome_conglomerate']['--th_finalAlignOverlap'] = args.th__ISinRefGenome_conglomerate__finalAlignOverlap[0]

		if self.dict_th['mergeLocalCounts']['--separator'] != self.dict_th['calcInRef_posISorient_v2']['--separator']: 
			sys.exit('Error: mergeLocalCounts.separator and calcInRef_posISorient_v2.separator should be the same.' )

		self.dict_th['appendEstimatedWrtRef']['--th_overlap'] = args.th__appendEstimatedWrtRef__overlapTh[0]


	def getThresholds_asList(self, scriptName):
		asList = [] 
		for key in self.dict_th[scriptName]: 
			asList.append(key) 
			asList.append(str(self.dict_th[scriptName][key]))

		# print (asList)
		# return self.dict_th[scriptName]
		return asList
		 

	def printThresholds(self): 
		print (self.dict_th)


#################################### TOP_LVL 
def runWaIS(dir_out, isRunSpades, fnList_assembly, fn_forwardReads, fn_reverseReads, fn_ISseqs, fn_reference, fn_referenceAnnotations, thresholds, keepTmp, path_to_script): 

	(dir_out_wais, dir_out_waisTmp, dir_out_waisFinal, dir_out_spades) = createOutputDirStruct(dir_out, isRunSpades)

	## Run spades, or set as provided assembly 
	fn_assembly = handleSpades(isRunSpades, fnList_assembly, dir_out_spades, fn_forwardReads, fn_reverseReads)
	logging.info('Assembly file finalized as ' + fn_assembly)
	
	## Pipline steps
	logging.info('Starting WaIS pipeline steps.')
	
	fn_blastdb_assembly = dir_out_waisTmp + 'contigs_blastdb'
	fn_blastRes_IStoContigs = dir_out_waisTmp + 'IStoContigs_blastRes'
	fn_reads1_fa = dir_out_waisTmp + extractReadName(fn_forwardReads) + '.fasta'
	fn_reads2_fa = dir_out_waisTmp + extractReadName(fn_reverseReads) + '.fasta'
	fn_blastdb_ISseqs = dir_out_waisTmp + 'ISseqs_blastdb'
	fn_1_to_IS_blastRes = dir_out_waisTmp + '1_to_IS_blastRes'
	fn_2_to_IS_blastRes = dir_out_waisTmp + '2_to_IS_blastRes'
	fn_flanks_1 = dir_out_waisTmp + 'flanks_1.fasta'
	fn_flanks_2 = dir_out_waisTmp + 'flanks_2.fasta' 
	fn_flanks_1_filtered = dir_out_waisTmp + 'flanks_1_filtered.fasta'
	fn_flanks_2_filtered = dir_out_waisTmp + 'flanks_2_filtered.fasta'
	fn_flanks1Filtered_to_contigs = dir_out_waisTmp +  'flanks_1_to_contigs_filtered'
	fn_flanks2Filtered_to_contigs = dir_out_waisTmp + 'flanks_2_to_contigs_filtered'
	fn_only1Ids_moved = dir_out_waisTmp + 'only1.txt'
	fn_only2Ids_moved = dir_out_waisTmp + 'only2.txt'
	fn_both_moved = dir_out_waisTmp + 'both.txt'
	fn_only1Seqs = dir_out_waisTmp + 'only1.fasta'
	fn_only2Seqs = dir_out_waisTmp + 'only2.fasta'
	fn_only1ToContigs_blastRes = dir_out_waisTmp + 'only1ToContig.blastRes'
	fn_only2ToContigs_blastRes = dir_out_waisTmp + 'only2ToContig.blastRes'
	fn_only1AndFlanks2 = dir_out_waisTmp + 'only1AndFlanks2'
	fn_only2AndFlanks1 = dir_out_waisTmp + 'only2AndFlanks1'
	fn_flanks_1_filtered_B = dir_out_waisTmp + 'flanks_1_filteredB' 
	fn_flanks_2_filtered_B = dir_out_waisTmp + 'flanks_2_filteredB'
	fn_flanks1FilteredB_to_contigs = dir_out_waisTmp + 'flanks_1_to_contigs_filteredB.blastRes'
	fn_flanks2FilteredB_to_contigs = dir_out_waisTmp + 'flanks_2_to_contigs_filteredB.blastRes'
	fn_ISinContigs_all = dir_out_waisFinal + 'ISinContigs_all.gff'
	fn_ISinContigs_merged = dir_out_waisFinal + 'ISinContigs_merged.gff'
	fn_ISinContigs_ignoreOrient = dir_out_waisFinal + 'ISinContigs_ignoreOrient.gff'
	fn_ISinContigs_ignoreIStype = dir_out_waisFinal + 'ISinContigs_ignoreIStype.gff'
	
	fn_contigEstimates_merged = dir_out_waisFinal + 'estimates_merged.txt'
	fn_contigEstimates_ignoreOrient = dir_out_waisFinal + 'estimates_ignoreOrient.txt'
	fn_contigEstimates_ignoreIStype = dir_out_waisFinal + 'estimates_ignoreIStype.txt'

	fn_estimates_singleRow_merged = dir_out_waisFinal + 'estimates_singleRow_merged.txt'
	fn_estimates_singleRow_ignoreOrient = dir_out_waisFinal + 'estimates_singleRow_ignoreOrient.txt'
	fn_estimates_singleRow_ignoreIStype = dir_out_waisFinal + 'estimates_singleRow_ignoreIStype.txt'

	fn_out_rmFlanks = dir_out_waisTmp + 'outfile_rmFlanks_whenOneDirFullAlign'
	fn_out_loadOnlyAndPntFasta = dir_out_waisTmp + 'outfile_loadOnlyAndPntFasta'
	fn_out_calcInContig_posISOrient = dir_out_waisTmp + 'outfile_calcInContig_posISOrient'
	fn_out_mergeLocalCounts = dir_out_waisTmp + 'outfile' 

	## Reference files:
	fn_blastdb_reference = dir_out_waisTmp + 'reference_blastdb'
	fn_contigsToRef_blastRes = dir_out_waisTmp + 'contigsToRef_blastRes'
	fn_contigsToRef_gff = dir_out_waisTmp + 'contigsToRef_gff'
	fn_out_convertBlastToGff = dir_out_waisTmp + 'outfile_convertBlastToGff'


	fn_ISinRef_blastRes = dir_out_waisTmp + 'ISinRef_blastRes'
	fn_ISinRef_gff = dir_out_waisFinal + 'ISinRef.gff'

	fn_IStoRef_gff_all = dir_out_waisFinal + 'IStoRef_all.gff'
	fn_IStoRef_gff_merged = dir_out_waisFinal + 'IStoRef_merged.gff'
	fn_IStoRef_gff_ignoreOrient = dir_out_waisFinal + 'IStoRef_ignoreOrient.gff'
	fn_IStoRef_gff_ignoreIStype = dir_out_waisFinal + 'IStoRef_ignoreIStype.gff'
	fn_out = dir_out_waisTmp + 'outfile'
	fn_presAbs_merged = dir_out_waisFinal + 'inRef_presenceAbsence_merged'
	fn_presAbs_ignoreOrient = dir_out_waisFinal + 'inRef_presenceAbsence_ignoreOrient'
	fn_presAbs_ignoreIStype = dir_out_waisFinal + 'inRef_presenceAbsence_ignoreIStype'
	
	fn_foundNotfound = dir_out_waisFinal + 'inRef_found_notFound'
	
	## Prokka 
	### 1. Contigs: 
	dir_prokka_assembly = dir_out_wais + 'PROKKA_assembly'
	fn_prokka_assembly = '' 

	### 2. Reference: 
	dir_prokka_reference = dir_out_wais + 'PROKKA_reference'
	fn_prokka_reference = '' 
	
	## 3. RefAnnots
	fn_refAnnotsFN = '' 
	if fn_referenceAnnotations != None: 
		fn_refAnnotsFN = os.path.basename(fn_referenceAnnotations)

		# print (fn_refAnnotsFN)
	
	fn_allInterrupAnnots_merged = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'all_interupAnnot_merged'
	fn_allInterrupAnnots_ignoreOrient = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'all_interupAnnot_ignoreOrient'
	fn_allInterrupAnnots_ignoreIStype = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'all_interupAnnot_ignoreIStype'

	fn_onlyInterupAnnots_merged = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'only_interupAnnot_merged'
	fn_onlyInterupAnnots_ignoreOrient = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'only_interupAnnot_ignoreOrient'
	fn_onlyInterupAnnots_ignoreIStype = dir_out_waisFinal + fn_refAnnotsFN + '-' + 'only_interupAnnot_ignoreIStype'
	
	# sys.exit(); 

	# 1. makeblastdb of assembly
	makeBlastDb(fn_assembly, fn_blastdb_assembly, 'nucl')
	
	# 2. blast ISseqs to assembly
	doBlastn_outfmt7(fn_ISseqs, fn_blastdb_assembly, fn_blastRes_IStoContigs)

	
	# 3. Uncompress reads
	runSeqTk(fn_forwardReads, fn_reads1_fa)
	runSeqTk(fn_reverseReads, fn_reads2_fa)
	
	# 4. Makeblastdb of ISseqs
	makeBlastDb(fn_ISseqs, fn_blastdb_ISseqs, 'nucl')
	
	# 5. Blast reads to ISseqs 
	doBlastn_outfmt7(fn_reads1_fa, fn_blastdb_ISseqs, fn_1_to_IS_blastRes)
	doBlastn_outfmt7(fn_reads2_fa, fn_blastdb_ISseqs, fn_2_to_IS_blastRes)
	
	# 6. Get flanking seqs. (i.e. segment of read that does not align to ISseq)
	getFlankingSeqs(path_to_script, fn_1_to_IS_blastRes, fn_reads1_fa, fn_flanks_1, thresholds)
	getFlankingSeqs(path_to_script, fn_2_to_IS_blastRes, fn_reads2_fa, fn_flanks_2, thresholds)
	
	# 7. Refine flanks - remove pairs with complete alignment
	rmFlanks_whenOneDirFullAlign_v2(path_to_script, fn_1_to_IS_blastRes, fn_2_to_IS_blastRes, fn_flanks_1, fn_flanks_2, fn_flanks_1_filtered, fn_flanks_2_filtered, fn_out_rmFlanks, dir_out_waisTmp, thresholds)

	
	# 8. Blast filtered_flanks to contigs 
	doBlastn_outfmt7(fn_flanks_1_filtered, fn_blastdb_assembly, fn_flanks1Filtered_to_contigs)
	doBlastn_outfmt7(fn_flanks_2_filtered, fn_blastdb_assembly, fn_flanks2Filtered_to_contigs)

	# 9. Extract sequences from the onlyIds.
	loadOnlyAndPntFasta(path_to_script, fn_only1Ids_moved, fn_only2Ids_moved, fn_reads1_fa, fn_reads2_fa, fn_only1Seqs, fn_only2Seqs, fn_out_loadOnlyAndPntFasta)

	# 10. Blast 
	doBlastn_outfmt7(fn_only1Seqs, fn_blastdb_assembly, fn_only1ToContigs_blastRes)
	doBlastn_outfmt7(fn_only2Seqs, fn_blastdb_assembly, fn_only2ToContigs_blastRes)


	# 11. summarizeOnlyToContig 
	summarizeOnlyToContig(path_to_script, fn_only1ToContigs_blastRes, fn_flanks2Filtered_to_contigs, fn_only1AndFlanks2, thresholds)
	summarizeOnlyToContig(path_to_script, fn_only2ToContigs_blastRes, fn_flanks1Filtered_to_contigs, fn_only2AndFlanks1, thresholds)


	# 12. rmFlanks_whenPairMismatchContig
	rmFlanks_whenPairMismatchContig(path_to_script, fn_only1AndFlanks2, fn_only2AndFlanks1, fn_flanks_1_filtered, fn_flanks_2_filtered, fn_flanks_1_filtered_B, fn_flanks_2_filtered_B) 

	
	# 13. Blast
	doBlastn_outfmt7(fn_flanks_1_filtered_B, fn_blastdb_assembly, fn_flanks1FilteredB_to_contigs)
	doBlastn_outfmt7(fn_flanks_2_filtered_B, fn_blastdb_assembly, fn_flanks2FilteredB_to_contigs)

	# 14. calcInContig_posISOrient
	calcInContig_posISOrient(path_to_script, fn_blastRes_IStoContigs, fn_flanks1FilteredB_to_contigs, fn_flanks2FilteredB_to_contigs,  fn_ISinContigs_all, fn_out_calcInContig_posISOrient, thresholds) 

	
	# 15a. In Contigs: mergeLocalCounts (merge-overlaps)
	mergeLocalCounts(path_to_script, fn_ISinContigs_all, fn_contigEstimates_merged, fn_ISinContigs_merged, fn_out_mergeLocalCounts, [], thresholds, dir_out, fn_estimates_singleRow_merged)

	# 15b. In Contigs: mergeLocalCounts (merge, ignoreOrient,)
	mergeLocalCounts(path_to_script, fn_ISinContigs_all, fn_contigEstimates_ignoreOrient, fn_ISinContigs_ignoreOrient, fn_out_mergeLocalCounts, ['--ignoreOrient', 'True'], thresholds, dir_out, fn_estimates_singleRow_ignoreOrient)

	# 15c. In Contigs: mergeLocalCounts (merge, ignoreOrient, ignoreIStype)
	mergeLocalCounts(path_to_script, fn_ISinContigs_all, fn_contigEstimates_ignoreIStype, fn_ISinContigs_ignoreIStype, fn_out_mergeLocalCounts, ['--ignoreOrient', 'True', '--ignoreIStype', 'True'], thresholds, dir_out, fn_estimates_singleRow_ignoreIStype)

	
	## 16. If reference 
	if fn_reference != None:

		# 1. Make reference blast db. 
		makeBlastDb(fn_reference, fn_blastdb_reference, 'nucl')

		# 2. Blast IS to refs, and convert blastRes to gff3 with different levels of merges (and write to estimates as a 'golden' count).
		doBlastn_outfmt7(fn_ISseqs, fn_blastdb_reference, fn_ISinRef_blastRes) 
		convertBlastToGff_IStoRef(path_to_script, fn_ISinRef_blastRes, fn_contigEstimates_merged, fn_ISinRef_gff, fn_out_convertBlastToGff, [], fn_estimates_singleRow_merged) # converted (with merged)
		convertBlastToGff_IStoRef(path_to_script, fn_ISinRef_blastRes, fn_contigEstimates_ignoreOrient, fn_ISinRef_gff,fn_out_convertBlastToGff, ['--ignoreOrient', 'True'], fn_estimates_singleRow_ignoreOrient) # converted (with ignoreOrient)
		convertBlastToGff_IStoRef(path_to_script, fn_ISinRef_blastRes, fn_contigEstimates_ignoreIStype, fn_ISinRef_gff, fn_out_convertBlastToGff, ['--ignoreIStype', 'True'], fn_estimates_singleRow_ignoreIStype) # converted (with ignoreIStype) 

		# 3. Blast contigs to ref.
		doBlastn_outfmt7(fn_assembly, fn_blastdb_reference, fn_contigsToRef_blastRes) # blast contigs to refs

		# 4. (not used)
		convertBlastToGff(path_to_script, fn_contigsToRef_blastRes, fn_contigsToRef_gff, fn_out_convertBlastToGff) 
		# addRefCountsToEstimatesFile() (Script ready) (write to the estimates file)

		# 5. Calculate IS to Contigs (using mappings from IStoContigsToRef_blastResults )
		# calcInRef_posISorient_v2(fn_blastRes_IStoContigs, fn_ISinRef_blastRes, fn_contigsToRef_blastRes, fn_ISinContigs_merged, fn_IStoRef_gff_all, fn_out, [], thresholds) 

		## TODO: Redirect outputs from the following 3 to /dev/null (and delete).
		calcInRef_posISorient_v2(path_to_script, fn_blastRes_IStoContigs, fn_ISinRef_blastRes, fn_contigsToRef_blastRes, fn_ISinContigs_merged, fn_IStoRef_gff_merged, fn_out, ['--isMerged', 'True'], thresholds) 
		calcInRef_posISorient_v2(path_to_script, fn_blastRes_IStoContigs, fn_ISinRef_blastRes, fn_contigsToRef_blastRes, fn_ISinContigs_ignoreOrient, fn_IStoRef_gff_ignoreOrient, fn_out, ['--isMerged', 'True', '--ignoreOrient', 'True'], thresholds) 
		calcInRef_posISorient_v2(path_to_script, fn_blastRes_IStoContigs, fn_ISinRef_blastRes, fn_contigsToRef_blastRes, fn_ISinContigs_ignoreIStype, fn_IStoRef_gff_ignoreIStype, fn_out, ['--isMerged', 'True', '--ignoreOrient', 'True', '--ignoreIStype', 'True'], thresholds) 

		# 6. Append IS in ref. to estimates.txt files. 
		appendEstimatedWrtRef(path_to_script, fn_IStoRef_gff_merged, fn_ISinContigs_merged, fn_contigEstimates_merged, fn_estimates_singleRow_merged, fn_out)
		appendEstimatedWrtRef(path_to_script, fn_IStoRef_gff_ignoreOrient, fn_ISinContigs_ignoreOrient, fn_contigEstimates_ignoreOrient, fn_estimates_singleRow_ignoreOrient, fn_out)
		appendEstimatedWrtRef(path_to_script, fn_IStoRef_gff_ignoreIStype, fn_ISinContigs_ignoreIStype, fn_contigEstimates_ignoreIStype, fn_estimates_singleRow_ignoreIStype, fn_out)
		
		# 7. Generate a presence absence table with regards to the ref-IS (additionally those inserted at the new positions in the ref-genome). 
		genPresAbsTblWrtRef(path_to_script, fn_ISinRef_blastRes, fn_IStoRef_gff_merged, fn_presAbs_merged, dir_out)
		genPresAbsTblWrtRef(path_to_script, fn_ISinRef_blastRes, fn_IStoRef_gff_ignoreOrient, fn_presAbs_ignoreOrient, dir_out)
		genPresAbsTblWrtRef(path_to_script, fn_ISinRef_blastRes, fn_IStoRef_gff_ignoreIStype, fn_presAbs_ignoreIStype, dir_out)


		# 7. If user provides ref annotations: determine if those interrupted or not. 
		if fn_referenceAnnotations != None: 
			insertionsWrtRefAnnotations(path_to_script, fn_reference, fn_referenceAnnotations, fn_IStoRef_gff_merged, fn_allInterrupAnnots_merged, fn_onlyInterupAnnots_merged) 
			
			insertionsWrtRefAnnotations(path_to_script, fn_reference, fn_referenceAnnotations, fn_IStoRef_gff_ignoreOrient, fn_allInterrupAnnots_ignoreOrient, fn_onlyInterupAnnots_ignoreOrient) 
		
			insertionsWrtRefAnnotations(path_to_script, fn_reference, fn_referenceAnnotations, fn_IStoRef_gff_ignoreIStype, fn_allInterrupAnnots_ignoreIStype, fn_onlyInterupAnnots_ignoreIStype) 
		
		# insertionsWrtRefAnnotations()

		# Estimates of IStoRef.py
		# insertionsWrtRefAnnotations.py
	
		# genPresAbsTblWrtRef(fn_IStoRef_blastRes, fn_IStoRef_gff_all, fn_presenceAbsence)
	
		# ISinRefGenome_conglomerate(fn_IStoRef_blastRes, fn_IStoRef_gff_all, fn_foundNotfound, thresholds)

	"""	
	## If prokka 
	if isProkka: 
		# For contigs 
		runProkka(fn_assembly, dir_prokka_assembly, 'assembly')
		# getInterruptedAnnotationIds()
		# bedtoolsIntersect()

	
		if fn_reference != None: 
			runProkka(fn_reference, dir_prokka_reference, 'reference') 
			# bedtoolsIntersect()
	# ... 
	"""
	if fn_referenceAnnotations != None: 
		# Run a ISinsertionWrtAnnotations.py script
		pass 

	if not keepTmp: 
		shutil.rmtree(dir_out_waisTmp)


	
#################################### AUX - Prokka 
def bedtoolsIntersect():
	pass 

def runProkka(fn_input, dir_output, str_runningOn):
	command = ["prokka", fn_input, "--outdir", dir_output]

	# subprocess.run(command, shell=True) 

	runTheCommand(command, 'Running PROKKA on ' +  str_runningOn + '.')

	fn_prokka_out = glob.glob(dir_output + '/*.gff') 
	
	# print ("The output prokka filename is " + str(fn_prokka_out))
	return fn_prokka_out[0]

#################################### AUX - Calling WaIS scripts
def insertionsWrtRefAnnotations(path_to_script, fn_reference, fn_refAnnotations, fn_ISannotations, fn_allAnnots, fn_onlyInterupAnnots):
	command = ["python3", path_to_script + "/scripts/insertionsWrtRefAnnotations.py", "--reference", fn_reference, "--ref_annotation", fn_refAnnotations, '--IS_annotation', fn_ISannotations, '--outfile', fn_allAnnots, '--outfile_onlyInterrupted', fn_onlyInterupAnnots]

	# command = command + thresholds.getThresholds_asList('ISinRefGenome_conglomerate')
	# command = command + params 

	runTheCommand(command, 'Mapping identified IS positions to reference genome.') 
	

def appendEstimatedWrtRef(path_to_script, fn_ISmappedToRef, fn_ISinContigs, fn_estimates, fn_estimates_singleRow, fn_out): 
	command = ["python3", path_to_script + "/scripts/appendEstimatedWrtRef_v2.py", "--IStoRef", fn_ISmappedToRef, "--fn_estimates", fn_estimates, '--fn_estimates_singleRow', fn_estimates_singleRow, '--ISinContigs', fn_ISinContigs]

	# command = command + thresholds.getThresholds_asList('appendEstimatedWrtRef')
	# command = command + params 

	runTheCommand_redirectOutputToFile(command, 'Mapping identified IS positions to reference genome.', fn_out) 


def ISinRefGenome_conglomerate(fn_IStoRef_blastRes, fn_IStoRef_gff, fn_res, thresholds):
	command = ["python3", "wais/scripts/ISinRefGenome_conglomerate.py", "--IStoRef_blastRes", fn_IStoRef_blastRes, "--IStoRef_gff", fn_IStoRef_gff] # + " > " + fn_res] 

	command = command + thresholds.getThresholds_asList('ISinRefGenome_conglomerate')

	runTheCommand_redirectOutputToFile(command, 'Getting IS in reference, conglomerated.', fn_res)

	
	# subprocess.run(command, shell=True)


def genPresAbsTblWrtRef(path_to_script, fn_ISinRef_blastRes, fn_IStoRef_gff, fn_presenceAbsence, isolateId):
	command = ["python3", path_to_script + "/scripts/genPresAbsTblWrtRef.py", "--IStoRef_blast", fn_ISinRef_blastRes, "--IStoRef_mapped", fn_IStoRef_gff, "--isolateId", isolateId] #  + " > " + fn_presenceAbsence]

	runTheCommand_redirectOutputToFile(command, 'Generating IS site presence, or absence, for an isolate w.r.t. reference.', fn_presenceAbsence)

	# subprocess.run(command, shell=True) 

def calcInRef_posISorient_v2(path_to_script, fn_blastRes_IStoContigs, fn_IStoRef_blastRes, fn_contigToRef_blastRes, fn_ISinContigs_all, fn_IStoRef_all_gff, fn_out, params, thresholds): 
	command = ['python3', path_to_script + '/scripts/calcInRef_posISorient_v2.py', '--IStoContig', fn_blastRes_IStoContigs, '--directIStoRef', fn_IStoRef_blastRes, '--contigToRef', fn_contigToRef_blastRes, '--IS_annotation', fn_ISinContigs_all, '--fn_out', fn_IStoRef_all_gff] # ' ' + params + ' > ' + fn_out]

	command = command + thresholds.getThresholds_asList('calcInRef_posISorient_v2')
	command = command + params

	runTheCommand_redirectOutputToFile(command, 'Mapping identified IS positions to reference genome.', fn_out) 
	# subprocess.run(command, shell=True)
	
def convertBlastToGff_IStoRef(path_to_script, fn_IStoRef_blastRes, fn_estimates, fn_IStoRef_gff, fn_out_convertBlastToGff, additionalArgs, fn_estimates_singleRow): 
	command = ['python3', path_to_script + '/scripts/convertBlastToGff_ref.py', '--blastRes', fn_IStoRef_blastRes, '--fn_estimates', fn_estimates, '--out', fn_IStoRef_gff, '--fn_estimates_singleRow', fn_estimates_singleRow]

	command = command + additionalArgs

	runTheCommand_redirectOutputToFile(command, 'Converting contigs-to-ref BLAST results to gff3.', fn_out_convertBlastToGff)
	 

def convertBlastToGff(path_to_script, fn_contigsToRef_blastRes, fn_contigsToRef_gff, fn_out_convertBlastToGff):
	command = ['python3', path_to_script + '/scripts/convertBlastToGff.py', '--blastRes', fn_contigsToRef_blastRes, '--out', fn_contigsToRef_gff] # , ' > ' + fn_out_convertBlastToGff]

	runTheCommand_redirectOutputToFile(command, 'Converting contigs-to-ref BLAST results to gff3.', fn_out_convertBlastToGff)

	# subprocess.run(command, shell=True)

def mergeLocalCounts(path_to_script, fn_ISinContig_all, fn_contigEstimates, fn_ISinContigs_merged, fn_out, additionalParams, thresholds, isolateId, fn_estimates_singleRow): 
	command = ['python3', path_to_script + '/scripts/mergeLocalCounts.py', '--fn_ISinGff', fn_ISinContig_all, '--fnOut_estimates', fn_contigEstimates,  '--fnOut_gff3_merged', fn_ISinContigs_merged, '--isolateId', isolateId, '--fnOut_estimates_singleRow', fn_estimates_singleRow] #  ' ' + additionalParams + ' > ' + fn_out] 

	command = command + thresholds.getThresholds_asList('mergeLocalCounts')
	command = command + additionalParams


	runTheCommand(command, 'Merging local counts - to get estimates of IS insertions.')

	# subprocess.run(command, shell=True)

def calcInContig_posISOrient(path_to_script ,fn_blastRes_IStoContigs, fn_flanks1FilteredB, fn_flanks2FilteredB, fn_ISinContig_gff, fn_out, thresholds): 
	command = ['python3', path_to_script + '/scripts/calcInContig_posISOrient.py', '--direct', fn_blastRes_IStoContigs, '--flanks1ToContigs', fn_flanks1FilteredB, '--flanks2ToContigs', fn_flanks2FilteredB, '--output_gff', fn_ISinContig_gff] #, ' > ' + fn_out]

	command = command + thresholds.getThresholds_asList('calcInContig_posISOrient')
	# subprocess.run(command, shell=True) 
	runTheCommand(command, 'Calculating IS insertions in contigs.')

def rmFlanks_whenPairMismatchContig(path_to_script, fn_only1AndFlanks2, fn_only2AndFlanks1, fn_flanks1Filtered, fn_flanks2Filtered, fn_flanks1FilteredB, fn_flanks2FilteredB):

	command = ['python3', path_to_script + '/scripts/rmFlanks_whenPairMismatchContig.py', '--only2AndFlanks1', fn_only2AndFlanks1, '--only1AndFlanks2', fn_only1AndFlanks2, '--flanks1', fn_flanks1Filtered, '--flanks2', fn_flanks2Filtered, '--out_flanks1', fn_flanks1FilteredB, '--out_flanks2', fn_flanks2FilteredB]

	runTheCommand(command, 'Removing pairs when mismatched IS reads.')


def summarizeOnlyToContig(path_to_script, fn_only1ToContigs_blastRes, fn_flanks2Filtered_to_contigs, fn_only1AndFlanks2, thresholds):
	command = ["python3", path_to_script + "/scripts/summarizeOnlyToContig.py", "--only", fn_only1ToContigs_blastRes, "--flanks", fn_flanks2Filtered_to_contigs] #, " > " + fn_only1AndFlanks2]

	command = command + thresholds.getThresholds_asList('summarizeOnlyToContig')
	
	runTheCommand_redirectOutputToFile(command, 'Summarising "only" blastRes to contigs.',fn_only1AndFlanks2)

	# subprocess.run(command, shell=True) 


def loadOnlyAndPntFasta(path_to_script, fn_only1, fn_only2, fn_reads1, fn_reads2, fn_out_reads1, fn_out_reads2, fn_out):

	command = ['python3', path_to_script + '/scripts/loadOnlyAndPntFasta.py', '--only1', fn_only1, '--only2', fn_only2, '--reads1', fn_reads1, '--reads2', fn_reads2, '--outfile_reads1', fn_out_reads1, '--outfile_reads2', fn_out_reads2] # , ' > ' + fn_out]

	runTheCommand(command, 'Getting the fasta sequences of the filtered (1) flanks.')
	# subprocess.run(command, shell=True) 


def rmFlanks_whenOneDirFullAlign_v2(path_to_script, blastRes_1, blastRes_2, flanks_1, flanks_2, fn_out_flanks1Filtered, fn_out_flanks2Filtered, fn_out, dir_out_waisTmp, thresholds):
	command1 = ["python3", path_to_script + "/scripts/rmFlanks_whenOneDirFullAlign_v2.py", "--blastRes_1", blastRes_1, "--blastRes_2", blastRes_2, "--flanks_1", flanks_1, "--flanks_2", flanks_2, "--out_flanks_1", fn_out_flanks1Filtered, "--out_flanks_2", fn_out_flanks2Filtered] #, " > " + fn_out]

	command1 = command1 + thresholds.getThresholds_asList('rmFlanks_whenOneDirFullAlign_v2')

	command2 = ["mv", "only1.txt", "only2.txt", "both.txt", dir_out_waisTmp]

	runTheCommand(command1, 'Filtering flanks part 1 (removing flanks, when pair is complete IS')
	runTheCommand(command2, 'Moving extracted flanks to tmp dir')

	# subprocess.run(command1, shell=True)
	# subprocess.run(command2, shell=True)


def getFlankingSeqs(path_to_script, reads_to_IS_blastRes, reads, fn_out, thresholds):
	command = ['python3', path_to_script + '/scripts/getFlankingSeqs.py', '--reads_to_IS', reads_to_IS_blastRes, '--reads', reads] #  + ' > ' + fn_out]
	
	command = command + thresholds.getThresholds_asList('getFlankingSeqs')
	
	runTheCommand_redirectOutputToFile(command, 'WaIS - extracting flanking sequences.', fn_out)
	# subprocess.run(command, shell=True)
#################################### AUX 
def extractReadName(fn_read): 

	fn_ = os.path.basename(fn_read)
	arr = fn_.split('.')

	return arr[0]


#################################### AUX - SEQtk
def runSeqTk(fn_reads, fn_out): 

	command = ['seqtk', 'seq', '-A', '-N', fn_reads, '>', fn_out]

	# subprocess.run(, shell=True)
	runTheCommand_redirectOutputToFile(command, 'Decompressing ' + fn_reads + ' to fasta using seqtk ' + fn_out, fn_out)

#################################### AUX - BLAST
def makeBlastDb(input, output, dbtype): 
	
	
	command = ["makeblastdb",  "-in", input, "-out", output, "-dbtype", dbtype]
	# subprocess.run("makeblastdb -in " + input + " -out " + output + " -dbtype " + dbtype, shell=True)
	runTheCommand(command, 'Making BLAST database ' + str(output))


def runTheCommand(command, printStr):
	logging.info('START: ' + printStr)
	logging.info('Command: ' + ' '.join(command))

	# subprocess.check_call(command)
	cmd = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	cmdOut, cmdErr = cmd.communicate()

	if cmdOut: 
		cmdOut = cmdOut.decode("utf-8").replace('\\n', '\n')
		logging.info(cmdOut) 

	if cmdErr or cmd.returncode != 0: 
		logging.info('ERROR: ' + printStr)
		logging.error(cmdErr)
		sys.stderr.write('ERROR: ' + printStr)
		sys.exit('See ' + str(logFile) + ' for details.')


	logging.info('COMPLETE: ' + printStr + '\n\n')

	
def runTheCommand_redirectOutputToFile(command, printStr, fn_out):

	fh_out = open(fn_out, 'w+')

	logging.info('START: ' + printStr)
	logging.info('Command: ' + ' '.join(command))

	# subprocess.check_call(command)
	cmd = subprocess.Popen(command, stdout=fh_out, stderr=subprocess.STDOUT)
	cmdOut, cmdErr = cmd.communicate()

	if cmdOut: 
		cmdOut = cmdOut.decode("utf-8").replace('\\n', '\n')
		logging.info(cmdOut) 

	if cmdErr or cmd.returncode != 0: 
		logging.info('ERROR: ' + printStr)
		logging.error(cmdErr)
		sys.stderr.write('ERROR: ' + printStr)
		sys.exit('See ' + str(logFile) + ' for details.')


	logging.info('COMPLETE: ' + printStr + '\n\n')

def doBlastn_outfmt7(query, db, output):

	command = ['blastn', '-task', 'megablast', '-query', query, '-db', db, '-outfmt', "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen", '-out', output]

	runTheCommand(command, 'BLASTing query ' + query + ' against db ' + db)

#################################### AUX - SPADES
def handleSpades(isRunSpades, fn_assembly, dir_out_spades, fn_reads1, fn_reads2): 

	if isRunSpades == False: 
		return fn_assembly[0]

	# Run spades and return assembly filename 

	logging.info('Running spades.py')

	fh_spades_shell_stdout = open(dir_out_spades + 'shell_stdout', 'w+')
	fh_spades_shell_stderr = open(dir_out_spades + 'shell_stderr', 'w+')

	subprocess.run('spades.py --tmp-dir /tmp --pe1-1 ' + fn_reads1 + ' --pe1-2 ' + fn_reads2 +  ' -o ' + dir_out_spades, stdout=fh_spades_shell_stdout, shell=True, stderr=fh_spades_shell_stderr)

	fh_spades_shell_stdout.close() 
	fh_spades_shell_stderr.close() 


	fn_spades_assembly = dir_out_spades + 'contigs.fasta' 

	if not os.path.exists(fn_spades_assembly):
		sys.exit("\nError: Spades did not generated a contigs.fasta file." + '\n\n')
	# sys.exit('\nTodo: Run spades and return assembly filename (+ spades arguments)' + '\n\n')
	return fn_spades_assembly

#################################### AUX - DIR
def createOutputDirStruct(dir_out, isRunSpades): 

	dir_out_wais = dir_out + 'WaIS/'
	dir_out_waisTmp = dir_out_wais + 'Tmp/'
	dir_out_waisFinal = dir_out_wais + 'Final/'
	dir_out_spades = dir_out + 'Spades/'

	os.mkdir(dir_out)
	os.mkdir(dir_out_wais)
	os.mkdir(dir_out_waisTmp)
	os.mkdir(dir_out_waisFinal)

	if isRunSpades == True: 
		os.mkdir(dir_out_spades)

	return (dir_out_wais, dir_out_waisTmp, dir_out_waisFinal, dir_out_spades)

#################################### MAIN 

def checkSpades(runSpades, assembly): 

	if ((runSpades == False) and assembly == None) or (runSpades == True and assembly): 
		sys.exit("\nError: exactly one of the following options is required:\n --runSpades\n --assembly <contigs.fasta>\n\n")


def checkIfFileExists(filename): 		
	if not os.path.exists(filename): 
		sys.exit('\nError: the specified assembly file does not exist ' + filename + '\n')

def returnSlashIfMissing(dirName): 

	if not re.match('/$', dirName): 
		return '/'
	
	return ''


def main(): 



	
	parser = argparse.ArgumentParser(description='Determine where the insertion sequences (IS) are in a reference genome, or an assembly (generated using short-reads) - using short-read sequences.')

	parser.add_argument('--outputDir', required=True, nargs=1)
	parser.add_argument('--runSpades', action='store_true', help='Assemble genome as part of WaIS using SPAdes.')
	parser.add_argument('--spadesOptions', nargs=1, default=[''], help="Options to send to spades.")
	parser.add_argument('--assembly', nargs=1, help='')
	parser.add_argument('--ISseqs', required=True, nargs=1, help='Fasta file containing the IS sequences to find.') 
	parser.add_argument('--reads_1', required=True, nargs=1, help="Illumina reads forward file.")
	parser.add_argument('--reads_2', required=True, nargs=1, help="Illumina reads reverse file.")
	
	parser.add_argument('--keepTmp', action='store_true', help='Keep the WaisTmp folder. Default=False.')
	
	## Reference 
	parser.add_argument('--reference', nargs=1, help='A reference genome to map the insertion sequences identified in the assembly.', default=[None])

	## Prokka 
	# parser.add_argument('--prokka', action='store_true', help='Run PROKKA to annotate genome and identify those distrupted by IS.')
	parser.add_argument('--referenceAnnotations', nargs=1, help='Annotations of genes, or other regions of interest, in the .gbk (genbank) format.', default=[None])

	## Spades options
	# --path_to_spades 
	# --options_to_spades

	## WaIS thresholds
	### getFlankingSeqs.py (Step 7.)
	parser.add_argument('--th_getFlankingSeqs_overlap', type=int, nargs=1, metavar='INT', help="Number of base pairs of IS alignment to include. Default=0.", default=[0])
	parser.add_argument('--th_getFlankingSeqs_minChoppedLen', type=int, nargs=1, metavar='INT', help="Minimun chopped sequence length (excluding overlapping sequence) to keep that sequence for further analysis. Default=18.", default=[18])

	### rmFlanks_whenOneDirFullAlign_v2 (Step 8.)
	parser.add_argument('--th__rmFlanks_whenOneDirFullAlign_v2__minChoppedLen', type=int, nargs=1, metavar='INT', help="Minimun chopped sequence length (excluding overlapping sequence) to keep that sequence for further analysis. Default=18.", default=[18])
	parser.add_argument('--th__rmFlanks_whenOneDirFullAlign_v2__buffer', type=int, nargs=1, metavar='INT', help="Buffer sequence length of the alignment of SKESA to complete genome to include, to see if the IS lies nearby. Default=0.", default=[0])

	### summarizeOnlyToContig.py (Step 11.)
	parser.add_argument('--th_summarizeOnlyToContig_minPident', type=int, nargs=1, metavar='INT', help="Minimum percent identity of alignment to keep (flanks to contig). Default=0.", default=[0])
	parser.add_argument('--th_summarizeOnlyToContig_alignAndLenDiff', type=int, nargs=1, metavar='INT', help="Alignment and length difference to keep (flanks to contig). Default=10000.", default=[10000])

	### calcInContig_posISOrient (Step 14.) 
	parser.add_argument('--th__calcInContig_posISOrient__minPident', type=int, nargs=1, metavar='INT', help='Minimum percent identity of alignment to keep (flanks to contig). Default=0.', default=[0]) 
	parser.add_argument('--th__calcInContig_posISOrient__minPalignLen', type=int, nargs=1, metavar='INT', help='The minimum percentage length of flank-sequence that aligns with contig. Default=0.', default=[0])
	parser.add_argument('--th__calcInContig_posISOrient__minAlignLen', type=int, nargs=1, metavar='INT', help='The minimum alignment length. Default=18.', default=[18]) 
	parser.add_argument('--th__calcInContig_posISOrient__minPident_direct', type=int, nargs=1, metavar='INT', help='Minimum percent identity of alignment to keep (direct IS to contig). Default=95.', default=[95]) 
	parser.add_argument('--th__calcInContig_posISOrient__minAlignLen_direct', type=int, nargs=1, metavar='INT', help='The minimum alignment length of direct IS to contig. Default=18.', default=[18]) 
	parser.add_argument('--th__calcInContig_posISOrient__clus_start', type=int, nargs=1, metavar='INT', help='Minimum number of clusters to evaluate in kMeans for each contig. Default=2.', default=[2])
	parser.add_argument('--th__calcInContig_posISOrient__clus_end', type=int, nargs=1, metavar='INT', help='Maximum number of clusters to evaluate in kMeans for each contig. Default=11.', default=[11]) 
	parser.add_argument('--th__calcInContig_posISOrient__minFlankDepth', type=int, nargs=1, metavar='INT', help='The minimum number of flanks (/reads) to keep IS position. Default=10.', default=[10]) 

	### mergeLocalCounts (Steps 15a, 15b, 15c)
	parser.add_argument('--th_mergeLocalCounts_toCountAsEdge', type=int, nargs=1, metavar='INT', help='This distance from the (start+th) or (end-th) is checked to determine if insertion is counted as being an edge insertion or an insertion in the middle. Default = 100 (bps).', default=[100]) 
	parser.add_argument('--th_mergeLocalCounts_forMergingOverlaps', type=int, nargs=1, metavar='INT', help='Distance between two IS, before calling them merges, Default = 20.', default=[20]) 
	parser.add_argument('--th_mergeLocalCounts_separator', nargs=1, metavar='STRING', help='Separate unresolvable IS and orientations. Default is ":". Choose a separator that is not an alphabet in the IS identifier.', default=[':']) 

	### (Reference:) calcInRef_posISorient_v2.py 
	parser.add_argument('--th__calcInRef_posISorient_v2__merge', type=int, nargs=1, metavar='INT', help='Threshold (i.e. sequence length in base pairs) to merge overlapping IS. Value only used in conjection with --isMerged. Default=20.', default=[20])
	parser.add_argument('--th__calcInRef_posISorient_v2__alignDiff_IStoContig', type=int, metavar='INT', help="Maximum alignment difference between align-length and contig-length. Default=15.", default=[5])
	parser.add_argument('--th__calcInRef_posISorient_v2__separator', nargs=1, metavar='STRING', help='Separate unresolvable IS and orientations. Default is ":". Choose a separator that is not an alphabet in the IS identifier.', default=[':']) 

	### (Reference:) ISinRefGenome_conglomerate.py
	parser.add_argument('--th__ISinRefGenome_conglomerate__toMergePosFound', type=int, nargs=1, metavar='INT', help='To merge IS-positions found within \'th\'. Default=20.', default=[20])
	parser.add_argument('--th__ISinRefGenome_conglomerate__finalAlignOverlap', type=int, nargs=1, metavar='INT', help='To merge IS-positions found within \'th\'. Default=0.', default=[0])

	parser.add_argument('--th__appendEstimatedWrtRef__overlapTh', type=int, nargs=1, metavar='INT', help='Distance in basepairs to describe as single IS insertion. Default=20.', default=[20])

	parser.usage = parser.format_help()
	args = parser.parse_args()
	
	thresholds = Thresholds(args)
	# thresholds.printThresholds()
	# print (thresholds.getThresholds_asList('getFlankingSeqs'))

	# print ("Printing the args:")
	# print (args)

	checkSpades(args.runSpades, args.assembly)
	
	## Checking if input/output files exist
	checkIfFileExists(args.ISseqs[0])
	checkIfFileExists(args.reads_1[0])
	checkIfFileExists(args.reads_2[0])
	if args.assembly and len(args.assembly) > 0: 
		checkIfFileExists(args.assembly[0])

	if args.reference[0] != None: 
		checkIfFileExists(args.reference[0])

	if args.referenceAnnotations[0] != None and args.reference[0] == None:
		sys.exit('\nError: you have provided a referenceAnnotations file, please remove this option or also provide a reference genome (option --reference).' + '\n') 
		# checkIfFileExists(args.reference[0])

	elif args.referenceAnnotations[0] != None and args.reference[0] != None: 
		checkIfFileExists(args.referenceAnnotations[0])

	if os.path.exists(args.outputDir[0]): 
		sys.exit('\nError: the folder already exists ' + args.outputDir[0] + '\n')

	## Adding slash
	
	args.outputDir[0] = args.outputDir[0] + returnSlashIfMissing(args.outputDir[0])	

	## Checking input thresholds

	path_to_script = re.split(r'\/', os.path.realpath(__file__))
	path_to_script = '/'.join(path_to_script[:-1])
	print(path_to_script)


	runWaIS(args.outputDir[0], args.runSpades, args.assembly, args.reads_1[0], args.reads_2[0], args.ISseqs[0], args.reference[0], args.referenceAnnotations[0], thresholds, args.keepTmp, path_to_script)




if __name__=='__main__':
	main() 