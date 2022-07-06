#!/usr/bin/python3
import argparse
import re
import sys

from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import silhouette_score

from scipy import stats
import matplotlib.pyplot as plt

from gap_statistic import OptimalK
from unidip import UniDip
import ckwrap 
from ckmeans import ckmeans

import contextlib
import os

import matplotlib.pyplot as plt
import random

import warnings
warnings.filterwarnings("error")

PROG_NAME = 'WiIS'
SO_ = "insertion_sequence" # Sequence ontology

####################################### TOP_LVL
def doTheCalculation(fn_direct, fn_flanks1ToContigs, fn_flanks2ToContigs, th_minPident, th_minPalignLen, th_minAlignLen, numClusStart, numClusEnd, th_minPident_direct, th_minAlignLen_direct, fn_out_gff, th_minFlankDepth):

    list_allContigsEncountered = []

    # 1. Load the direct IS_to_contigs blast results
    dict_direct = loadDirectIStoContig(fn_direct, list_allContigsEncountered, th_minPident_direct, th_minAlignLen_direct)




    # 1. Load flanks_to_contigs blast results
    dict_flanksToContigs = dict() # dict_{(contigId, contigLen)} => [(flankId, ISid, orient_flankToIS, flankToContig_start, flankToContig_end, isContigFlipped, orient_flankToContig, alignLen)]

    loadFlanks(dict_flanksToContigs, fn_flanks1ToContigs, th_minPident, th_minPalignLen, th_minAlignLen)
    loadFlanks(dict_flanksToContigs, fn_flanks1ToContigs, th_minPident, th_minPalignLen, th_minAlignLen)



    # 2. Extract positions for each contig, and cluster.
    # with open(fn_out_vcf, 'w+') as fh_out:
    
    dict_colors = {} # dict_[ISid] => "color" 

    fh_gff3 = open(fn_out_gff, 'w+')
    
    # print ('##gff-version 3.1.26')

    contigsEncountered = [] 
    for contigInfo in dict_flanksToContigs:

        
        # Determine the positions (using kmeans, silhoutte index and ckmeans)
        (dict_eachFlankInsPos, dict_mergedClusters) =  clusterSplitAndMerge_posCalc(dict_flanksToContigs[contigInfo], numClusStart, numClusEnd)


        if contigInfo in dict_direct:

            # Merge any flankInsPos to direct
            mergeToDirect(dict_direct[contigInfo], dict_eachFlankInsPos)

            if len(dict_mergedClusters) > 0: 
                mergeToDirect(dict_direct[contigInfo], dict_mergedClusters)
                
            # Noting that it was added 
            addTo_contigsEncountered(contigsEncountered, contigInfo)



                
            # Output
            for (ISid, alignPosInContig_start, alignPosInContig_end, orientISToContig, isContigFlipped) in dict_direct[contigInfo]: 

                if ISid not in dict_colors: 
                    dict_colors[ISid] = generateRandomColor()

                printGff3Line(fh_gff3, contigInfo[0], PROG_NAME, SO_, alignPosInContig_start, alignPosInContig_end, "Inf", orientISToContig, '.', "Name=" + ISid + ';' + 'color=' + dict_colors[ISid])
                # print (contigInfo[0] + "\t" + "Direct" + "\t "+ str(aDirectIS))

        

        if len(dict_mergedClusters) > 0: 
            for clusNum in dict_mergedClusters:
                for (ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect) in dict_mergedClusters[clusNum]:

                    if (flankCounts) > th_minFlankDepth and idx_mergeWithDirect == -1 : 
                        
                        # Noting that it was added 
                        addTo_contigsEncountered(contigsEncountered, contigInfo)
                
                        # Output 
                        doTheFlankOutput(fh_gff3, contigInfo[0], ISid, dict_orient, dict_side, flankCounts, alignLens, positions, dict_colors)
                        

        else: 
            for clusNum in dict_eachFlankInsPos:
                for (ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect) in dict_eachFlankInsPos[clusNum]:
                    if (flankCounts) > th_minFlankDepth and idx_mergeWithDirect == -1 : 

                        # Noting that it was added 
                        addTo_contigsEncountered(contigsEncountered, contigInfo)
                
                        # Output 
                        doTheFlankOutput(fh_gff3, contigInfo[0], ISid, dict_orient, dict_side, flankCounts, alignLens, positions, dict_colors)

    
    for contigInfo in dict_direct: 
        if contigInfo not in contigsEncountered: 
            for (ISid, alignPosInContig_start, alignPosInContig_end, orientISToContig, isContigFlipped) in dict_direct[contigInfo]: 

                # Noting that it was added 
                addTo_contigsEncountered(contigsEncountered, contigInfo)
                
                if ISid not in dict_colors: 
                    dict_colors[ISid] = generateRandomColor()

                printGff3Line(fh_gff3, contigInfo[0], PROG_NAME, SO_, alignPosInContig_start, alignPosInContig_end, "Inf", orientISToContig, '.', "Name=" + ISid + ';' + 'color=' + dict_colors[ISid])



    fh_gff3.close() 


    # Reload file and add Gff3 header (gff3 version + sequence_regions)
    with open(fn_out_gff, 'r') as original: 
        data = original.read()
    
    with open(fn_out_gff, 'w') as modified: 
        modified.write('##gff-version 3.1.26' + '\n')
        for (contigId, contigLen) in contigsEncountered: 
            modified.write('##sequence-region ' + contigId + ' ' + '1' + ' ' + str(contigLen) + '\n')
        
        modified.write(data)
    # Print contig lengths in gff file header.

    """
    print();
    if len(dict_mergedClusters) > 0: 
        print('\nMerged 2 (merged clusters):')
        for clusNum in dict_mergedClusters:
            for (ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect) in dict_mergedClusters[clusNum]:

                if (flankCounts) > th_minAlignLen: 
                    avgAlignLen = np.average(np.asarray(alignLens))
                    minPos = np.amin(positions)
                    maxPos = np.amax(positions)
                    modePos = stats.mode(positions)

                    print (str(clusNum) + "\t" + str(ISid) + "\t" + str(dict_orient) + "\t" + str(dict_side) + "\t" + str(flankCounts) + "\t" + str(minPos) + ":" + str(maxPos) + "\t" + str(modePos) + "\t" + str(idx_mergeWithDirect))
    
    else: 
        for clusNum in dict_eachFlankInsPos:
            for (ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect) in dict_eachFlankInsPos[clusNum]:
                if (flankCounts) > th_minAlignLen: 
                    avgAlignLen = np.average(np.asarray(alignLens))
                    minPos = np.amin(positions)
                    maxPos = np.amax(positions)
                    modePos = stats.mode(positions)

                    print (str(clusNum) + "\t" + str(ISid) + "\t" + str(dict_orient) + "\t" + str(dict_side) + "\t" + str(flankCounts) + "\t" + str(minPos) + ":" + str(maxPos) + "\t" + str(modePos) + "\t" + str(idx_mergeWithDirect))

    print ('\n\n\n')
    # break;
    """

def generateRandomColor(): 
    color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    print(color)
    
    return color 


def printGff3Line(fh_, refSeq, source, featureType, startPos, endPos, score, strand, phase, attributes):
    fh_.write (refSeq + '\t' + source + '\t' + featureType + "\t" + str(startPos) + "\t" + str(endPos) + "\t" + str(score) + "\t" + strand + "\t" + phase + "\t" + attributes + '\n')
    



def addTo_contigsEncountered(contigsEncountered, contigInfo): 

    if contigInfo not in contigsEncountered: 
        contigsEncountered.append(contigInfo) 




def clusterSplitAndMerge_posCalc(flanksToAContig, numClusStart, numClusEnd):


    # C-means clustering
    (k_centers, labels, silhouette_coefficients) = calculateClusters(flanksToAContig, numClusStart, numClusEnd)

    # print(silhouette_coefficients)
    # Calculate best cluster
    # print ("The best cluster is " + str(silhouette_coefficients.index(max(silhouette_coefficients)) + numClusStart))
    
    if len(silhouette_coefficients) == 0:
        ## print ("Make all labels 1")
        dict_flankCounts = clus_posAndInfo_split(np.ones(len(flanksToAContig)), flanksToAContig)
        dict_eachFlankInsPos = mergeSplitsWithinClus(dict_flankCounts)
    
    else:
        idx_bestCluster = silhouette_coefficients.index(max(silhouette_coefficients)) # The best cluster number #k - numClusStart
        ## print ('Best cluster num.: ' + str(idx_bestCluster + numClusStart) + "; Silhouette coeff.:" + str(silhouette_coefficients[idx_bestCluster]))
        # print(k_centers[idx])
        dict_flankCounts = clus_posAndInfo_split(labels[idx_bestCluster], flanksToAContig)

        # print (dict_flankCounts)
        dict_eachFlankInsPos = mergeSplitsWithinClus(dict_flankCounts) # When position is same.

    # print (dataset[:,0])
    
    (bestK_ckmeans, uniqFreq) = segmentIn1D(dict_eachFlankInsPos, 1, numClusEnd)

    bestK_kmeans = len(dict_flankCounts)
    inferredClusNum = 0 
    for key in dict_flankCounts: 
        inferredClusNum = inferredClusNum + len(dict_flankCounts[key])

    is1Dbetween = False

    if bestK_ckmeans >= bestK_kmeans and bestK_ckmeans <= inferredClusNum: 
        is1Dbetween = True

    ## print('Ckmeans: ' + str(bestK_ckmeans) + "; UniqFreq: " + str(uniqFreq) + "; ClusNum: " + str(bestK_kmeans) + "; InferredClusNum: " + str(inferredClusNum) + "; is1Dbetween: " + str(is1Dbetween))

    dict_mergedClusters = {}
    if is1Dbetween == False and uniqFreq < bestK_kmeans: 
        # Time for a second merge
        ## print ('Time for a second merge')
        dict_mergedClusters = merge2BasedOnISid(dict_eachFlankInsPos)


    return (dict_eachFlankInsPos, dict_mergedClusters)
    

def segmentIn1D(dict_eachFlankInsPos, numClusStart, numClusEnd):

    # Gather all the (1D) positions in the contig
    allContigPos = []

    for clusNum in dict_eachFlankInsPos:
        for (ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect) in dict_eachFlankInsPos[clusNum]:
            allContigPos = allContigPos + positions



    # Minimum is 1, max is num. of uni elements; if num. of uniq elemenst is less than user max.

    allContigPos_npArr = np.asarray(allContigPos)
    uniqFreq = len(np.unique(allContigPos_npArr))

    if numClusEnd > uniqFreq: 
        numClusEnd = uniqFreq
    
    # print (allContigPos)

    #res = ckwrap.cksegs(allContigPos_npArr) # (numClusStart, numClusEnd))
    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        res = ckmeans(allContigPos_npArr, (numClusStart, numClusEnd))
    

    """ GapTest
    dataset = x.astype(np.float)

    optimalK = OptimalK(parallel_backend=None, clusterer=special_clustering_func)
    n_clusters = optimalK(dataset, cluster_array=np.arange(1, 10))

    print ('Gap statistic: ' + str(n_clusters))
    """ 
    """ Dip test
    dataset = np.msort(allContigPos)
    intervals = UniDip(allContigPos).run()

    print ('Dip-test-intervals: ' + str(intervals))
    """

    return (res.k, uniqFreq)


def special_clustering_func(X, k):
    """
    Special clustering function which uses the MeanShift
    model from sklearn.

    These user defined functions *must* take the X and a k
    and can take an arbitrary number of other kwargs, which can
    be pass with `clusterer_kwargs` when initializing OptimalK
    """

    # Here you can do whatever clustering algorithm you heart desires,
    # but we'll do a simple wrap of the MeanShift model in sklearn.


    kmeans_kwargs = {
        "init": "random",
        "n_init": 10,
        "max_iter": 300
    }

    kmeans = KMeans(n_clusters=k, **kmeans_kwargs).fit(X)


    #m = MeanShift()
    #m.fit(X)

    # Return the location of each cluster center,
    # and the labels for each point.
    return kmeans.cluster_centers_, kmeans.predict(X)

################################################# OUTPUT
def doTheFlankOutput(fh_, contigId, ISid, dict_orient, dict_side, flankCounts, alignLens, positions, dict_colors): 
    avgAlignLen = np.average(np.asarray(alignLens))
    medianAlignLens = np.median(alignLens)
    meanAlignLens = np.mean(alignLens) 
    minPos = np.amin(positions)
    maxPos = np.amax(positions)
    modePos = stats.mode(positions)
    medianPos = np.median(positions)
    meanPos = np.mean(positions) 



    gff3_startPos = minPos
    gff3_endPos = maxPos 
    gff3_attributes = 'Name=' + ISid + ';' 
    gff3_orient = '.' 

    if len(modePos[0]) > 0:
        gff3_startPos = modePos[0][0] 
        gff3_endPos = modePos[0][0]
        gff3_attributes = gff3_attributes + 'modeCount=' + str(modePos[1][0]) + ";"

    if (maxPos != minPos):     
        gff3_attributes = gff3_attributes + 'meanPos=' + "{:.1f}".format(meanPos) + ";"
        gff3_attributes = gff3_attributes + 'medianPos=' + str(medianPos) + ";"
    gff3_attributes = gff3_attributes + 'maxPos=' + str(maxPos) + ";"
    gff3_attributes = gff3_attributes + 'minPos=' + str(minPos) + ";"
    

    if len(dict_orient) == 1: 
        for key in dict_orient: 
            gff3_orient = key
    else: 
        for key in dict_orient: 
            gff3_attributes = gff3_attributes + 'orient_' + key + '=' + str(dict_orient[key]) + ';'

    for key in dict_side:
        gff3_attributes = gff3_attributes + "side_" + key + '=' + str(dict_side[key]) + ';'

    gff3_attributes = gff3_attributes + "meanAlignLen=" + "{:.1f}".format(meanAlignLens) + ';'

    gff3_attributes = gff3_attributes + "medianAlignLen=" + "{:.1f}".format(medianAlignLens) + ';'
    
    if ISid not in dict_colors: 
        dict_colors[ISid] = generateRandomColor()
    
    gff3_attributes = gff3_attributes + 'color=' + dict_colors[ISid]

    printGff3Line (fh_, contigId, PROG_NAME, SO_, gff3_startPos, gff3_endPos, flankCounts, gff3_orient, '.', gff3_attributes)

def printToVcf(fh_out, flankOrDirect, contigId, ISid, positions, orient_IStoContig, flankDepth, avgFlankDepth):
    #fh_out.write(contigId + "\t" + str(positions[0]) + ":" + str(positions[1]))

    #fh_out.write('\n')
    pass

#################################################
def genTheStatsForAPos(arrPositions):
    npArr = np.asarray(arrPositions)

    minVal = np.amin(npArr)
    maxVal = np.amax(npArr)
    rangeMinMax = np.ptp(npArr)
    medianVal = np.median(npArr)
    modeVal = stats.mode(npArr)
    percentile = np.percentile(npArr, 50)

    ## print ('\t\t' + str(minVal) + "\t" + str(maxVal) + "\t" +  str(medianVal) + "\t" + str(modeVal) + "\t" + str(rangeMinMax))
    ## print('\t\tPercentile_50:' + str(percentile))

def calcSideOfIS(flankId, orient_flankToContig):
    Col_prePost = 1; Col_orient_IStoFlank = 7;

    arr_flankInfo = flankId.split("__")

    sideOfIS = None
    if orient_flankToContig == '+':
        if arr_flankInfo[Col_prePost] == 'Pre':
            sideOfIS = 'l'
        elif arr_flankInfo[Col_prePost] == 'Post':
            sideOfIS = 's'

    if orient_flankToContig == '-':
        if arr_flankInfo[Col_prePost] == 'Pre':
            sideOfIS = 's'
        elif arr_flankInfo[Col_prePost] == 'Post':
            sideOfIS = 'l'

    return sideOfIS

def mergeToDirect(directContigInfo, dict_eachFlankInsPos):


    idx_direct = 0
    for (ISid, alignPosInContig_start, alignPosInContig_end, orient_IStoContig, isContigFlipped) in directContigInfo:

        # print (str(alignPosInContig_start) + "\t" + str(alignPosInContig_end))
        for clusNum in dict_eachFlankInsPos:
            idx_flankInsPos = 0
            for (ISid_, dict_orient_, dict_side_, flankCounts_, alignLens_, positions_, idx_mergeWithDirect) in dict_eachFlankInsPos[clusNum]:

                flanksMin = np.amin(positions_)
                flanksMax = np.amax(positions_)

                if ISid == ISid_ and ((flanksMin >= (alignPosInContig_start - 1) and flanksMin <= (alignPosInContig_end + 1)) or (flanksMax >= (alignPosInContig_start - 1) and flanksMax <= (alignPosInContig_end + 1)) or (alignPosInContig_start >= (flanksMin - 1) and alignPosInContig_start <= (flanksMax + 1)) or (alignPosInContig_end >= (flanksMin - 1) and alignPosInContig_end <= (flanksMax + 1))):

                    idx_mergeWithDirect = idx_direct
                    dict_eachFlankInsPos[clusNum][idx_flankInsPos] = (ISid_, dict_orient_, dict_side_, flankCounts_, alignLens_, positions_, idx_mergeWithDirect)


                idx_flankInsPos = idx_flankInsPos + 1

        idx_direct = idx_direct + 1

def mergeSplitsWithinClus(dict_flankCounts):

    dict_eachInsPos = {} # dict_eachInsPos{clusNum} => [(ISid, dict_orient_IStoContig{orient}=>cnt, dict_sideOfIS{side}=>cnt, flankCnt, [alignLen], [insertPos1, insertPos2, ...], idx_mergeWithDirect)]
    # dataType = mode, median, mean, min, max, range, percentile25, percentile50, percentile100

    for clusNum in dict_flankCounts:
        for (ISid, orient_IStoContig, sideIS) in dict_flankCounts[clusNum]:

            # Check and add to dict
            if sideIS == 's':
                addToDict_eachInsPos(dict_eachInsPos, clusNum, ISid, orient_IStoContig, sideIS, dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][0], dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][1], dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][2])
            elif sideIS == 'l':
                addToDict_eachInsPos(dict_eachInsPos, clusNum, ISid, orient_IStoContig, sideIS, dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][0], dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][1], dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][3])

            ## print(str(clusNum) + "\t" + ISid + "\t" + orient_IStoContig + "\t" + sideIS + "\t" + str(dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][0]) + "\t" + str(min(dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][2])) + ":" + str(max(dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][2])) + "\t" + str(min(dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][3])) + ":" + str(max(dict_flankCounts[clusNum][(ISid, orient_IStoContig, sideIS)][3])))
            #else:
                # Check if needs to be merged

    return dict_eachInsPos

def merge2BasedOnISid(dict_old):
    dict_new = dict() # dict_{clusNum} => [(IS, ISid, dict_orient, dict_side, flankCounts, alignLens, positions, idx_mergeWithDirect), ...] 

    # Perhaps easier to do a manual merge (?)
    for old_clusNum in dict_old: 
        for (old_ISid, old_dict_orient, old_dict_side, old_flankCounts, old_alignLens, old_positions, old_idx_mergeWithDirect) in dict_old[old_clusNum]:  
            addToDict_merge2(dict_new, old_clusNum, old_ISid, old_dict_orient, old_dict_side, old_flankCounts, old_alignLens, old_positions, old_idx_mergeWithDirect)

    return dict_new
             
def addToDict_merge2(dict_new, old_clusNum,  old_ISid, old_dict_orient, old_dict_side, old_flankCounts, old_alignLens, old_positions, old_idx_mergeWithDirect): 


    # if dictionary is empty, add as-is. 
    if len(dict_new.keys()) == 0: 
        dict_new[old_clusNum] = [(old_ISid, old_dict_orient, old_dict_side, old_flankCounts, old_alignLens, old_positions, old_idx_mergeWithDirect)]
    
    # else, check first.
    else: 
        isMerged = False
        for new_clusNum in dict_new: 
            idx_ = 0
            for (new_ISid, new_dict_orient, new_dict_side, new_flankCounts, new_alignLens, new_positions, new_idx_mergeWithDirect) in dict_new[new_clusNum]: 
                new_minVal = np.amin(new_positions)
                new_maxVal = np.amax(new_positions)

                old_minVal = np.amin(old_positions)
                old_maxVal = np.amax(old_positions)


                if (new_ISid == old_ISid) and ((new_minVal >= old_minVal and new_minVal <= old_maxVal) or (new_maxVal >= old_minVal and new_maxVal <= old_maxVal) or (old_minVal >= new_minVal and old_minVal <= new_maxVal) or (old_maxVal >= new_minVal and old_maxVal <= new_maxVal)): # Do the merge (add counts to prev.)
                    
                    # loop through each orient.
                    for old_orient in old_dict_orient: 
                        if old_orient not in new_dict_orient: 
                            new_dict_orient[old_orient] = 0

                        new_dict_orient[old_orient] = new_dict_orient[old_orient] + old_dict_orient[old_orient]

                    # loop through each side.
                    for old_side in old_dict_side: 
                        if old_side not in new_dict_side: 
                            new_dict_side[old_side] = 0

                        new_dict_side[old_side] = new_dict_side[old_side] + old_dict_side[old_side]


                    # Add flank counts
                    new_flankCounts = new_flankCounts + old_flankCounts

                    # append alignLens
                    new_alignLens = new_alignLens + old_alignLens

                    # append positions 
                    new_positions = new_positions + old_positions

                    dict_new[new_clusNum][idx_] = (new_ISid, new_dict_orient, new_dict_side, new_flankCounts, new_alignLens, new_positions, new_idx_mergeWithDirect)


                    isMerged = True
                idx_ = idx_ + 1

        if isMerged == False: 

            if old_clusNum not in dict_new: 
                dict_new[old_clusNum] = [] 
            
            dict_new[old_clusNum].append((old_ISid, old_dict_orient, old_dict_side, old_flankCounts, old_alignLens, old_positions, old_idx_mergeWithDirect))

def addToDict_eachInsPos(dict_eachInsPos, clusNum, ISid, orient, side, flankCounts, alignLens, positions):

    # dict_eachInsPos[clusNum].append()

    if clusNum not in dict_eachInsPos:
        dict_orient = {}
        dict_orient[orient] = flankCounts

        dict_side = {}
        dict_side[side] = flankCounts

        dict_eachInsPos[clusNum] = [(ISid, dict_orient, dict_side, flankCounts, alignLens, positions, -1)]

        return
    # Check and add

    isMerged = False
    idx_ = 0
    for (ISid_, dict_orient_, dict_side_, flankCounts_, alignLens_, positions_, idx_mergeWithDirect) in dict_eachInsPos[clusNum]:
        minVal_ = np.amin(positions_)
        maxVal_ = np.amax(positions_)

        minVal = np.amin(positions)
        maxVal = np.amax(positions)


        if (ISid == ISid_) and ((minVal_ >= minVal and minVal_ <= maxVal) or (maxVal_ >= minVal and maxVal_ <= maxVal) or (minVal >= minVal_ and minVal <= maxVal_) or (maxVal >= minVal_ and maxVal <= maxVal_)): # Do the merge (add counts to prev.)
            if orient not in dict_orient_:
                dict_orient_[orient] = flankCounts
            else:
                dict_orient_[orient] = dict_orient_[orient] + flankCounts

            if side not in dict_side_:
                dict_side_[side] = flankCounts
            else:
                dict_side_[side] = dict_side_[side] + flankCounts

            flankCounts_ = flankCounts_ + flankCounts

            alignLens_ = alignLens_ + alignLens
            positions_ = positions_ + positions

            dict_eachInsPos[clusNum][idx_] = (ISid_, dict_orient_, dict_side_, flankCounts_, alignLens_, positions_, -1)

            isMerged = True


        idx_ = idx_ + 1

    if isMerged == False:
        dict_orient = {}
        dict_orient[orient] = flankCounts

        dict_side = {}
        dict_side[side] = flankCounts

        dict_eachInsPos[clusNum].append((ISid, dict_orient, dict_side, flankCounts, alignLens, positions, -1))



def clus_posAndInfo_split(labels, arr_flanksToContigs):
    dict_flankCounts = {} # dict_[clusLabel] => {(ISid, orient_IStoContig, sideIS)} => [flankCnt, [alignLen1, alignLen2, ...], [contigStart1, ...], [contigEnd1, ...]]


    for i in range(0, len(arr_flanksToContigs)):


        (flankId, ISid, orient_flankToIS, flankToContig_start, flankToContig_end, isContigFlipped, orient_flankToContig, alignLen) = arr_flanksToContigs[i]


        # Calculate orient_IStoContig
        orient_IStoContig = None
        if (orient_flankToIS == '-' and orient_flankToContig == '-') or (orient_flankToIS == '+' and orient_flankToContig == '+'):
            orient_IStoContig = '+'
        elif (orient_flankToIS == '-' and orient_flankToContig == '+') or (orient_flankToIS == '+' and orient_flankToContig == '-'):
            orient_IStoContig = '-'

        # Calculate side of IS
        sideOfIS = calcSideOfIS(flankId, orient_flankToContig)



        # Add to dict
        if labels[i] not in dict_flankCounts:
            dict_flankCounts[labels[i]] = {}

        key = (ISid, orient_IStoContig, sideOfIS)
        if key not in dict_flankCounts[labels[i]]:
            dict_flankCounts[labels[i]][key] = [0, [], [], []]

        dict_flankCounts[labels[i]][key][0] =  dict_flankCounts[labels[i]][key][0] + 1

        dict_flankCounts[labels[i]][key][1].append(alignLen)


        dict_flankCounts[labels[i]][key][2].append(int(flankToContig_start))
        dict_flankCounts[labels[i]][key][3].append(int(flankToContig_end))
        # print (flankId)
        """
        for clusLabel in labels:
            # print(clusLabel)

            if clusLabel not in dict_flankCounts:
                dict_flankCounts[clusLabel] = 0

            dict_flankCounts[clusLabel] = dict_flankCounts[clusLabel] + 1
        """
    return dict_flankCounts

####################################### AUX - direct IS to contigs
def loadDirectIStoContig(fn_direct, list_allContigsEncountered, th_minPident_direct, th_minAlignLen_direct):

    dict_direct = {} # dict_{(contigId, contigLen)} => [(ISid, alignPosInContig_start, alignPosInContig_end, orientIStoContig, isContigFlipped)]

    Col_ISid = 0
    Col_contigId = 1
    Col_contigLen = 13

    Col_pIdent = 2
    Col_alignLen = 3

    Col_contigStart = 8
    Col_contigEnd = 9
    Col_ISstart = 6
    Col_ISend = 7

    with open (fn_direct, 'r') as fh:
        for line in fh:
            if not re.match('^\#', line):
                # print (line)

                line = line.strip()
                arr = line.split("\t")

                if ((float(arr[Col_pIdent]) - th_minPident_direct) >= 0.01 and int(arr[Col_alignLen]) >= th_minAlignLen_direct):

                    if (arr[Col_contigId] not in list_allContigsEncountered):
                        list_allContigsEncountered.append(arr[Col_contigId])

                    # sys.stdout.write (arr[Col_ISid] + "\t" + arr[Col_contigId] + "\t" + arr[Col_pIdent] + "\t" + arr[Col_alignLen] + "\t" + arr[Col_contigStart] + ":" + arr[Col_contigEnd] + '\t' + arr[Col_ISstart] + ":" + arr[Col_ISend])

                    contigStart = int(arr[Col_contigStart])
                    contigEnd = int(arr[Col_contigEnd])

                    ISstart = int(arr[Col_ISstart])
                    ISend = int(arr[Col_ISend])

                    orient = ''; isContigFlipped = False
                    if (contigStart > contigEnd and ISstart > ISend) or (contigStart < contigEnd and ISstart < ISend): # same direction
                        orient = '+'
                        # print('\t+')
                    else:
                        orient = '-'
                        # print('\t-')

                    if contigStart > contigEnd:
                        tmp = contigStart
                        contigStart = contigEnd
                        contigEnd = tmp

                        isContigFlipped = True


                    if (arr[Col_contigId], int(arr[Col_contigLen])) not in dict_direct:
                        dict_direct[(arr[Col_contigId], int(arr[Col_contigLen]))] = []

                    dict_direct[(arr[Col_contigId], int(arr[Col_contigLen]))].append((arr[Col_ISid], contigStart, contigEnd, orient, isContigFlipped))

    """
    for key in dict_direct:
        print (key)
        for val in dict_direct[key]:
            print ("\t" + str(val))
    """

    return dict_direct

####################################### AUX - positions for each contig & cluster
def calculateClusters(aligningFlanks, numClusStart, numClusEnd):


    dataset = None
    isFirst = True

    for (flankId, ISid, orient_flankToIS, flankToContig_start, flankToContig_end, isContigFlipped, orient_flankToContig, alignLen) in aligningFlanks:

        newRow = np.array([float(flankToContig_start), float(flankToContig_end)])
        if (isFirst == True):
            dataset = [newRow]
            isFirst = False
        else:
            dataset = np.append(dataset, [newRow], axis=0)

    kmeans_kwargs = {
        "init": "random",
        "n_init": 10,
        "max_iter": 300
    }


    # print (dataset)
    #dat = dataset[:,0]
    #dat = np.msort(dat)
    #intervals = UniDip(dat).run()
    #print(intervals)


    kmeans = KMeans(n_clusters=1, **kmeans_kwargs).fit(dataset)
    #print(kmeans.inertia_)

    #optimalK = OptimalK(parallel_backend='None')
    #n_clusters = optimalK(dataset, cluster_array=np.arange(1, 5))

    #print ('N cluster is ')
    #print (n_clusters)
    #print (optimalK.gap_df.head())

    # Run multiple clusters
    k_centers = []
    silhouette_coefficients = []
    labels = []
    for k in range(numClusStart, numClusEnd):
        try:
            kmeans = KMeans(n_clusters=k, **kmeans_kwargs).fit(dataset)
            #print(kmeans.inertia_)
        except:
            ## print("Catching bug here... k=" + str(k));
            break
        # print (kmeans.labels_)
        # print("Num clusters = " + str(k))
        # print (kmeans.cluster_centers_)
        k_centers.append(kmeans.cluster_centers_)
        score = silhouette_score(dataset, kmeans.labels_)
        silhouette_coefficients.append(score)
        labels.append(kmeans.labels_)

        # getNumReads_ISid_dir(kmeans.cluster_centers_, )

    return (k_centers, labels, silhouette_coefficients)


####################################### AUX
def loadFlanks(dict_flanksToContigs, fn_flanks1ToContigs, th_minPident, th_minPalignLen, th_minAlignLen):

    Col_flankId = 0
    Col_contigId = 1
    Col_pIdent = 2
    Col_alignLen = 3
    Col_qSeqLen = 12
    Col_contig_alignStart = 8
    Col_contig_alignEnd = 9
    Col_flank_alignStart = 6
    Col_flank_alignEnd = 7
    Col_contigLen = 13

    dict_flankCounts = {} # dict_{flankId} => count # count == number of times passed

    with open(fn_flanks1ToContigs, 'r') as fh:
        for line in fh:

            if (not re.match("^\#", line)):
                arr = line.split("\t")

                pAlignLen = (int(arr[Col_alignLen])/int(arr[Col_qSeqLen])) * 100

                if float(arr[Col_pIdent]) >= th_minPident and int(arr[Col_alignLen]) >= th_minAlignLen and pAlignLen >= th_minPalignLen:
                    # print (arr[Col_contigId] + "\t" + arr[Col_pIdent] + "\t" + arr[Col_alignLen] + "\t" + arr[Col_qSeqLen] + "\t" + 'Pass' + "\t" + str(pAlignLen) + "\t" + arr[Col_contig_alignStart] +":"+arr[Col_contig_alignEnd])

                    orient_flankToContig = calcAlignOrient(int(arr[Col_contig_alignStart]), int(arr[Col_contig_alignEnd]), int(arr[Col_flank_alignStart]), int(arr[Col_flank_alignEnd]))

                    addToDict_flanksToContigs(dict_flanksToContigs, arr[Col_contigId], arr[Col_flankId], int(arr[Col_contig_alignStart]), int(arr[Col_contig_alignEnd]), orient_flankToContig, int(arr[Col_contigLen]), int(arr[Col_alignLen]))

                    if arr[Col_flankId] not in dict_flankCounts:
                        dict_flankCounts[arr[Col_flankId]] = 0;

                    dict_flankCounts[arr[Col_flankId]] = dict_flankCounts[arr[Col_flankId]] + 1
                # else:
                #    print (arr[Col_contigId] + "\t" + arr[Col_pIdent] + "\t" + arr[Col_alignLen] + "\t" + arr[Col_qSeqLen] + "\t" + 'Remove' + "\t" + str(pAlignLen))


    """
    for flankId, count in sorted(dict_flankCounts.items(), key=lambda item: item[1]):
        print (flankId + "\t" + str(count))
    """


def addToDict_flanksToContigs(dict_flanksToContigs, contigId, flankId, alignPosInContig_start, alignPosInContig_end, orient_flankToContig, contigLen, alignLen):

    # Extract ISid, orient_flankToIS

    (ISid, orient_flankToIS) = getFlankToIS_ISandOrient(flankId)

    # Make start always smaller (for later clustering?)
    isContigFlipped = False
    if (alignPosInContig_start > alignPosInContig_end):
        tmp = alignPosInContig_start
        alignPosInContig_start = alignPosInContig_end
        alignPosInContig_end = tmp
        isContigFlipped = True


    if (contigId, contigLen) not in dict_flanksToContigs:
        dict_flanksToContigs[(contigId, contigLen)] = []

    dict_flanksToContigs[(contigId, contigLen)].append((flankId, ISid, orient_flankToIS, alignPosInContig_start, alignPosInContig_end, isContigFlipped, orient_flankToContig, alignLen))


def getFlankToIS_ISandOrient(flankId):
    Col_ISid = 2
    Col_orient = 7

    arr = flankId.split("__")

    return (arr[Col_ISid], arr[Col_orient])

def calcAlignOrient(contig_alignStart, contig_alignEnd, flank_alignStart, flank_alignEnd):

    if (contig_alignStart >= contig_alignEnd and flank_alignStart >= flank_alignEnd) or (contig_alignStart <= contig_alignEnd and flank_alignStart <= flank_alignEnd):
        # print (str(contig_alignStart) + ":" + str(contig_alignEnd) + "\t" + str(flank_alignStart) + ":" + str(flank_alignEnd) + "\t" + '+')
        return '+'
    else:
        # print (str(contig_alignStart) + ":" + str(contig_alignEnd) + "\t" + str(flank_alignStart) + ":" + str(flank_alignEnd) + "\t" + '-')
        return '-'



####################################### MAIN

def main():
    parser = argparse.ArgumentParser(description='Get IS positions, which IS, and orientation for each contig in VCF format.')

    parser.add_argument('--direct', nargs=1, required=True, help="Direct IS to contig blast alignments")
    parser.add_argument('--flanks1ToContigs', nargs=1, required=True, help='Flanks1 to contigs')
    parser.add_argument('--flanks2ToContigs', nargs=1, required=True, help='Flanks2 to contigs')

    ## Thresholds
    parser.add_argument('--th_minPident', nargs=1, default=[0], type=float, help='Minimum percent identity of alignment to keep (flanks to contig)')
    parser.add_argument('--th_minPalignLen', nargs=1, default=[0], type=float, help='The minimum percentage length of flank-sequence that aligns with contig.')
    parser.add_argument('--th_minAlignLen', nargs=1, default=[18], type=float, help='The minimum alignment length.')

    # For direct IS to contigs: Want to be a bit more generous with alignLen, but more strict with pIdent.
    parser.add_argument('--th_minPident_direct', nargs=1, default=[95], type=float, help='Minimum percent identity of alignment to keep (direct IS to contig)')
    parser.add_argument('--th_minAlignLen_direct', nargs=1, default=[18], type=int, help='The minimum alignment length of direct IS to contig.')


    parser.add_argument('--kmeans_clus_start', nargs=1, default=[2], type=int, help='Minimum number of clusters to evaluate in kMeans for each contig.')
    parser.add_argument('--kmeans_clus_end', nargs=1, default=[11], type=int, help='Maximum number of clusters to evaluate in kMeans for each contig.')



    # parser.add_argument()
    parser.add_argument('--output_gff', nargs=1, required=True, help='Output filename.vcf to output the IS structural variants to w.r.t. contigs.')

    parser.add_argument('--th_minFlankDepth', nargs=1, default=[10], type=int, help='The minimum number of flanks (/reads) to keep IS position.')

    args = parser.parse_args()

    doTheCalculation(args.direct[0], args.flanks1ToContigs[0], args.flanks2ToContigs[0], args.th_minPident[0], args.th_minPalignLen[0], args.th_minAlignLen[0], args.kmeans_clus_start[0], args.kmeans_clus_end[0], float(args.th_minPident_direct[0]), int(args.th_minAlignLen_direct[0]), args.output_gff[0], int(args.th_minFlankDepth[0]))


if __name__ == '__main__':
    main()
