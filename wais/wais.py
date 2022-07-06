#!/usr/bin/python3 

import argparse
import sys
import os 
import re
import subprocess 
import logging 
from datetime import datetime

logging.basicConfig(filename='wais_' + str(datetime.now()) + '.log', level=logging.DEBUG)

#################################### TOP_LVL 
def runWaIS(dir_out, isRunSpades, fnList_assembly, fn_forwardReads, fn_reverseReads, fn_ISseqs, fn_reference): 

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
	
	fn_contigEstimates_merged = 'estimates_contigs_merged.txt'
	fn_contigEstimates_ignoreOrient = 'estimates_contigs_ignoreOrient.txt'
	fn_contigEstimates_ignoreIStype = 'estimates_contigs_ignoreIStype.txt'

	fn_out_rmFlanks = dir_out_waisTmp + 'outfile_rmFlanks_whenOneDirFullAlign'
	fn_out_loadOnlyAndPntFasta = dir_out_waisTmp + 'outfile_loadOnlyAndPntFasta'
	fn_out_calcInContig_posISOrient = dir_out_waisTmp + 'outfile_calcInContig_posISOrient'
	fn_out_mergeLocalCounts = dir_out_waisTmp + 'outfile' 

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
	getFlankingSeqs(fn_1_to_IS_blastRes, fn_reads1_fa, fn_flanks_1)
	getFlankingSeqs(fn_2_to_IS_blastRes, fn_reads2_fa, fn_flanks_2)
	
	# 7. Refine flanks - remove pairs with complete alignment
	rmFlanks_whenOneDirFullAlign_v2(fn_1_to_IS_blastRes, fn_2_to_IS_blastRes, fn_flanks_1, fn_flanks_2, fn_flanks_1_filtered, fn_flanks_2_filtered, fn_out_rmFlanks, dir_out_waisTmp)

	# 8. Blast filtered_flanks to contigs 
	doBlastn_outfmt7(fn_flanks_1_filtered, fn_blastdb_assembly, fn_flanks1Filtered_to_contigs)
	doBlastn_outfmt7(fn_flanks_2_filtered, fn_blastdb_assembly, fn_flanks2Filtered_to_contigs)

	# 9. Extract sequences from the onlyIds.
	loadOnlyAndPntFasta(fn_only1Ids_moved, fn_only2Ids_moved, fn_reads1_fa, fn_reads2_fa, fn_only1Seqs, fn_only2Seqs, fn_out_loadOnlyAndPntFasta)

	# 10. Blast 
	doBlastn_outfmt7(fn_only1Seqs, fn_blastdb_assembly, fn_only1ToContigs_blastRes)
	doBlastn_outfmt7(fn_only2Seqs, fn_blastdb_assembly, fn_only2ToContigs_blastRes)

	# 11. summarizeOnlyToContig 
	summarizeOnlyToContig(fn_only1ToContigs_blastRes, fn_flanks2Filtered_to_contigs, str(95), str(40), fn_only1AndFlanks2)
	summarizeOnlyToContig(fn_only2ToContigs_blastRes, fn_flanks1Filtered_to_contigs, str(95), str(40), fn_only2AndFlanks1)

	# 12. rmFlanks_whenPairMismatchContig
	rmFlanks_whenPairMismatchContig(fn_only1AndFlanks2, fn_only2AndFlanks1, fn_flanks_1_filtered, fn_flanks_2_filtered, fn_flanks_1_filtered_B, fn_flanks_2_filtered_B) 

	# 13. Blast
	doBlastn_outfmt7(fn_flanks_1_filtered_B, fn_blastdb_assembly, fn_flanks1FilteredB_to_contigs)
	doBlastn_outfmt7(fn_flanks_2_filtered_B, fn_blastdb_assembly, fn_flanks2FilteredB_to_contigs)

	# 14. calcInContig_posISOrient
	calcInContig_posISOrient(fn_blastRes_IStoContigs, fn_flanks1Filtered_to_contigs, fn_flanks2Filtered_to_contigs, str(85), str(90), str(18), fn_ISinContigs_all, fn_out_calcInContig_posISOrient) 


	# 15a. In Contigs: mergeLocalCounts (merge-overlaps)
	mergeLocalCounts(fn_ISinContigs_all, fn_contigEstimates_merged, fn_ISinContigs_merged, str(20), fn_out_mergeLocalCounts, '')

	# 15b. In Contigs: mergeLocalCounts (merge, ignoreOrient,)
	mergeLocalCounts(fn_ISinContigs_all, fn_contigEstimates_merged, fn_ISinContigs_merged, str(20), fn_out_mergeLocalCounts, ' --ignoreOrient True ')

	# 15c. In Contigs: mergeLocalCounts (merge, ignoreOrient, ignoreIStype)
	mergeLocalCounts(fn_ISinContigs_all, fn_contigEstimates_ignoreOrient, fn_ISinContigs_ignoreOrient, str(20), fn_out_mergeLocalCounts, ' --ignoreOrient True --ignoreIStype True ')

	## If reference 
	# 16. calcInRef_posISorient_v2

	## If prokka 
	# 17. prokka 

	# ... 


#################################### AUX - Calling WaIS scripts
def mergeLocalCounts(fn_ISinContig_all, fn_contigEstimates, fn_ISinContigs_merged, th_forMergingOverlaps, fn_out, additionalParams): 
	command = 'python3 wais/scripts/mergeLocalCounts.py  --fn_ISinGff ' + fn_ISinContig_all + ' --fnOut_estimates ' + fn_contigEstimates + ' --th_forMergingOverlaps  ' + th_forMergingOverlaps + ' --fnOut_gff3_merged ' + fn_ISinContigs_merged + ' ' + additionalParams + ' > ' + fn_out 

	subprocess.run(command, shell=True)

def calcInContig_posISOrient(fn_blastRes_IStoContigs, fn_flanks1FilteredB, fn_flanks2FilteredB, th_minPident, th_minPalignLen, th_minAlignLen, fn_ISinContig_gff, fn_out): 
	command = 'python3 wais/scripts/calcInContig_posISOrient.py --direct ' + fn_blastRes_IStoContigs + ' --flanks1ToContigs ' + fn_flanks1FilteredB + ' --flanks2ToContigs ' + fn_flanks2FilteredB + ' --th_minPident ' + th_minPident + ' --th_minPalignLen ' + th_minPalignLen + ' --th_minAlignLen ' + th_minAlignLen + ' --output_gff ' + fn_ISinContig_gff + ' > ' + fn_out

	subprocess.run(command, shell=True) 


def rmFlanks_whenPairMismatchContig(fn_only1AndFlanks2, fn_only2AndFlanks1, fn_flanks1Filtered, fn_flanks2Filtered, fn_flanks1FilteredB, fn_flanks2FilteredB):
	command = 'python3 wais/scripts/rmFlanks_whenPairMismatchContig.py --only2AndFlanks1 ' + fn_only2AndFlanks1 + ' --only1AndFlanks2 ' + fn_only1AndFlanks2 + ' --flanks1 ' + fn_flanks1Filtered + ' --flanks2 ' + fn_flanks2Filtered + ' --out_flanks1 ' + fn_flanks1FilteredB + ' --out_flanks2 ' + fn_flanks2FilteredB

	subprocess.run(command, shell=True)


def summarizeOnlyToContig(fn_only1ToContigs_blastRes, fn_flanks2Filtered_to_contigs, th_minPident, th_alignAndLenDiff, fn_only1AndFlanks2):
	command = "python3 /srv/scratch/lanlab/Sandeep/Wiis/Scripts/summarizeOnlyToContig.py --only " + fn_only1ToContigs_blastRes + " --flanks " + fn_flanks2Filtered_to_contigs + " --th_minPident " + th_minPident + " --th_alignAndLenDiff " + th_alignAndLenDiff + " > " + fn_only1AndFlanks2
	
	subprocess.run(command, shell=True) 


def loadOnlyAndPntFasta(fn_only1, fn_only2, fn_reads1, fn_reads2, fn_out_reads1, fn_out_reads2, fn_out):

	command = 'python3 wais/scripts/loadOnlyAndPntFasta.py --only1 ' + fn_only1 + ' --only2 ' + fn_only2 + ' --reads1 ' + fn_reads1 + ' --reads2 ' + fn_reads2 + ' --outfile_reads1 ' + fn_out_reads1 + ' --outfile_reads2 ' + fn_out_reads2 + ' > ' + fn_out

	subprocess.run(command, shell=True) 


def rmFlanks_whenOneDirFullAlign_v2(blastRes_1, blastRes_2, flanks_1, flanks_2, fn_out_flanks1Filtered, fn_out_flanks2Filtered, fn_out, dir_out_waisTmp):
	command1 = "python3 wais/scripts/rmFlanks_whenOneDirFullAlign_v2.py --blastRes_1 " + blastRes_1 + " --blastRes_2 " + blastRes_2 + " --flanks_1 " + flanks_1 + " --flanks_2 " + flanks_2 + " --out_flanks_1 " + fn_out_flanks1Filtered + " --out_flanks_2 " + fn_out_flanks2Filtered + " > " + fn_out
	command2 = "mv only1.txt only2.txt both.txt " + dir_out_waisTmp

	subprocess.run(command1, shell=True)
	subprocess.run(command2, shell=True)


def getFlankingSeqs(reads_to_IS_blastRes, reads, fn_out):
	command = 'python3 wais/scripts/getFlankingSeqs.py --reads_to_IS ' + reads_to_IS_blastRes + ' --reads ' + reads + ' > ' + fn_out
	
	subprocess.run(command, shell=True)
#################################### AUX 
def extractReadName(fn_read): 

	fn_ = os.path.basename(fn_read)
	arr = fn_.split('.')

	return arr[0]


#################################### AUX - SEQtk
def runSeqTk(fn_reads, fn_out): 
	
	logging.info('Decompressing ' + fn_reads + ' to fasta using seqtk ' + fn_out)

	subprocess.run('seqtk seq -A -N ' + fn_reads + ' > ' + fn_out, shell=True)

#################################### AUX - BLAST
def makeBlastDb(input, output, dbtype): 
	
	logging.info('Creating blastdb ' + output)

	subprocess.run("makeblastdb -in " + input + " -out " + output + " -dbtype " + dbtype, shell=True)
		
def doBlastn_outfmt7(query, db, output):

	logging.info('Blasting query ' + query + ' against db ' + db)

	subprocess.run('blastn -task megablast -query ' + query + ' -db ' + db + ' -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out ' + output, shell=True) 
 

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
def main(): 
	parser = argparse.ArgumentParser(description='Determine where the insertion sequences (IS) are in a reference genome, or an assembly (generated using short-reads) - using short-read sequences.')

	parser.add_argument('--outputDir', required=True, nargs=1)
	parser.add_argument('--runSpades', action='store_true', help='')
	parser.add_argument('--assembly', nargs=1, help='')
	parser.add_argument('--ISseqs', required=True, nargs=1, help='Fasta file containing the IS sequences to find.') 
	parser.add_argument('--reads_1', required=True, nargs=1, help="Illumina reads forward file.")
	parser.add_argument('--reads_2', required=True, nargs=1, help="Illumina reads reverse file.")
	parser.add_argument('--reference', nargs=1, help='A reference genome to map the insertion sequences identified in the assembly.', default=[None])
	# '--keepTmp'
	## Reference 


	## Spades options
	# --path_to_spades 
	# --options_to_spades

	## WaIS thresholds

	args = parser.parse_args()
	
	checkSpades(args.runSpades, args.assembly)
	
	## Checking if input/output files exist
	checkIfFileExists(args.ISseqs[0])
	checkIfFileExists(args.reads_1[0])
	checkIfFileExists(args.reads_2[0])
	if args.assembly and len(args.assembly) > 0: 
		checkIfFileExists(args.assembly[0])

	if os.path.exists(args.outputDir[0]): 
		sys.exit('\nError: the folder already exists ' + args.outputDir[0] + '\n')

	## Adding slash
	
	args.outputDir[0] = args.outputDir[0] + returnSlashIfMissing(args.outputDir[0])	

	## Checking input thresholds

	runWaIS(args.outputDir[0], args.runSpades, args.assembly, args.reads_1[0], args.reads_2[0], args.ISseqs[0], args.reference[0])




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

if __name__=='__main__':
	main() 