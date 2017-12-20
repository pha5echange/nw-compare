# nw_compare_a11.py
# Version a11
# by jmg - j.gagen*AT*gold*DOT*ac*DOT*uk
# Nov 16th 2017

# Licence: http://creativecommons.org/licenses/by-nc-sa/3.0/

# Makes inclusion lists from master lists and temporally-appropriate GEXFs
# Generates master dictionary of EN/MB-to-Wiki genre labels from 2015 lists
# Generates master dictionary of MB-to-Wiki genre labels from 2015 lists
# Utilises omegaYear input to find GEXFs
# Uses these to generate echo_inclusion_lists and mb_inclusion_lists
# Then uses dictionary to generate wiki_inclusion_lists


# encoding=utf8
import sys
reload(sys)
sys.setdefaultencoding('utf8')

# Import packages
import os
import resource
import numpy as np
import networkx as nx
import community
import matplotlib.pyplot as plt
from networkx.algorithms.approximation import clique
from collections import OrderedDict
from datetime import datetime

versionNumber = ("a11")

# Initiate timing of run
runDate = datetime.now()
startTime = datetime.now()

# Create 'logs' subdirectory if necessary
if not os.path.exists("logs"):
    os.makedirs("logs")

# Create 'results' subdirectories if necessary
if not os.path.exists("results"):
	os.makedirs("results")

if not os.path.exists("results/analysis"):
	os.makedirs("results/analysis")

if not os.path.exists("results/communities"):
	os.makedirs("results/communities")

# Create 'networks' subdirectory (for images) if necessary
if not os.path.exists("networks"):
    os.makedirs("networks")

# Begin
print ('\n' + "Inclusion List Maker | Version " + versionNumber + " | Starting...")

# Get user input
print
try:
	dateIP = int(input ("Enter an omega year for this run (`0' or `2015' if unsure): "))
except:
	pass

if dateIP == 0:
	dateIP = 2015

omegaYear = str(dateIP)

# Open file for writing log
logPath = os.path.join("logs", omegaYear + '_nw_compare_' + versionNumber  + '_' + str(runDate) + '_log.txt')
runLog = open(logPath, 'a')


runLog.write ("==========================================================================" + '\n' + '\n')
runLog.write ("Inclusion List Maker | Version " + versionNumber + '\n' + '\n')

# Open files for inclusion_list ouput
echoListPath = os.path.join("inclusion-lists", omegaYear + '_echo_inclusion_list.txt')
echoListFile = open(echoListPath, 'w')

mbListPath = os.path.join("inclusion-lists", omegaYear + '_mb_inclusion_list.txt')
mbListFile = open(mbListPath, 'w')

wikiListPath = os.path.join("inclusion-lists", omegaYear + '_wiki_inclusion_list.txt')
wikiListFile = open(wikiListPath, 'w')

# Open files for writing gexfs
wikiGexfPath = os.path.join("gexf", omegaYear + '_nw_compare_wiki_' + versionNumber + '.gexf')
wikiGexfFile = open(wikiGexfPath, 'w')

echoGexfPath = os.path.join("gexf", omegaYear + '_nw_compare_echo_' + versionNumber + '.gexf')
echoGexfFile = open(echoGexfPath, 'w')

mbGexfPath = os.path.join("gexf", omegaYear + '_nw_compare_mb_' + versionNumber + '.gexf')
mbGexfFile = open(mbGexfPath, 'w')

# Open file for `communities' ouput files
wikiCommPath = os.path.join("results/communities", omegaYear + '_nw_compare_wiki_' + versionNumber + '_' + str(runDate) + '_communities.txt')
wikiCommFile = open(wikiCommPath, 'w')

echoCommPath = os.path.join("results/communities", omegaYear + '_nw_compare_echo_' + versionNumber + '_' + str(runDate) + '_communities.txt')
echoCommFile = open(echoCommPath, 'w')

mbCommPath = os.path.join("results/communities", omegaYear + '_nw_compare_mb_' + versionNumber + '_' + str(runDate) + '_communities.txt')
mbCommFile = open(mbCommPath, 'w')

# Open file for analysis results
anPath = os.path.join("results/analysis", omegaYear + '_nw_compare_' + versionNumber  + '_' + str(runDate) + '_analysis.txt')
anFile = open(anPath, 'w')

# Open files to write images
wikiImgPath = os.path.join("networks", omegaYear + '_nw_compare_wiki_' + versionNumber + '_nw.eps')
#wikiImg = open (wikiImgPath, 'w')

echoImgPath = os.path.join("networks", omegaYear + '_nw_compare_echo_' + versionNumber + '_nw.eps')
#echoImg = open (echoImgPath, 'w')

mbImgPath = os.path.join("networks", omegaYear + '_nw_compare_mb_' + versionNumber + '_nw.eps')
#mbImg = open (mbImgPath, 'w')

anFile.write ("==========================================================================" + '\n' + '\n')
anFile.write ("Inclusion List Maker | Version " + versionNumber + '\n' + '\n')

# Uncomment the line below to facilitate optional self-loop removal
#selfLoopIP = int(input ("Enter 1 here to remove self-loop edges: "))
# Auto-remove self-loops - comment out the line below if self-loop removal is optional
selfLoopIP = 1

# Create dictionaries for en/mb-to-wiki and mb-to-wiki genre mapping
matchedNodesPath = os.path.join("inclusion-lists", 'en_wd_mb_matched.txt')
matchedNodesFile = open(matchedNodesPath, 'r')

print
print("Storing mapped nodes (Echo:Wiki and MB:Wiki) in Dictionaries...")
print

enwdMapDict = {}
mbwdMapDict = {}

for line in matchedNodesFile:
	echoNode, wikiNode, mbNode = line.split(",")
	echoNode = echoNode.strip().replace('\n',"")
	wikiNode = wikiNode.strip().replace('\n',"")
	mbNode = mbNode.strip().replace('\n',"")
	enwdMapDict[echoNode] = wikiNode
	mbwdMapDict [mbNode] = wikiNode
	print(str(echoNode) + ',' + str(wikiNode) + '\n')
	print(str(mbNode) + ',' + str(wikiNode) + '\n')

print("Dictionaries complete: ")
print("Echo:Wiki Dictionary: ")
print(enwdMapDict)
print("MB:Wiki Dictionary")
print(mbwdMapDict)

# Read the EchoNest GEXF and generate graph
print
print ("Reading echoNest GEXF file and generating graph... ")
enInputPath = os.path.join("gexf", omegaYear + '_echonest.gexf')
echoNW = nx.read_gexf(enInputPath)

# Calculate basic Echo graph statistics
print
print ('Calculating various things for Echo graph...' + '\n')
echoNodes = nx.number_of_nodes(echoNW)
echoEdges = nx.number_of_edges(echoNW)
echoDensity = nx.density(echoNW)
echoNodeList = nx.nodes(echoNW)
echoNodeList.sort()
echoSelfLoopTotal = echoNW.number_of_selfloops()
echoConnections = echoEdges - echoSelfLoopTotal

print ("Omega Year: " + omegaYear)
print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Self-loops: ' + str(echoSelfLoopTotal))
print ('Echo Connections (edges minus self-loops): ' + str(echoConnections))
print ('Echo Density: ' + str(echoDensity))
print
print (str(nx.info(echoNW)))
print
print (str(echoNodeList))

anFile.write ('Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Self-loops: ' + str(echoSelfLoopTotal) + '\n')
anFile.write ('Echo Connections (edges minus self-loops): ' + str(echoConnections) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n' + '\n')
anFile.write (str(nx.info(echoNW)) + '\n' + '\n')

runLog.write ('Echo nodes: ' + '\n' + '\n')
runLog.write (str(echoNodeList) + '\n' + '\n')

# Remove nodes that are not on inclusion-list (`inclusion-lists/echo_inclusion_list.txt')

# Open echo_inclusion_list
echoRemovedCount = 0
print ('\n' + 'Checking list and removing echo nodes...' + '\n')
runLog.write ('Checking list and removing echo nodes...' + '\n' + '\n')
anFile.write ('Checking list and removing echo nodes...' + '\n' + '\n')
	
echoListPath = os.path.join("inclusion-lists", 'master_echo_inclusion_list.txt')
echoList = open(echoListPath, 'r')

# Check for genres and remove if not present
echoIncludeNodes = [omegaYear + '_' + line.strip() for line in echoList]

for i in echoNodeList:
	if not i in echoIncludeNodes:
		echoNW.remove_node(i) 
		print ('Removed echo node ' + str(i))
		runLog.write ('Removed echo node ' + str(i) + '\n')
		# anFile.write (str(i) + '\n')
		echoRemovedCount += 1
	else:
		echoNodeOP = str(i).replace(omegaYear + "_","")
		echoListFile.write(echoNodeOP + '\n')
		wikiNodeOP = enwdMapDict[echoNodeOP]
		wikiListFile.write(wikiNodeOP + '\n')

wikiListFile.close()

print ('\n'+ "Removed " + str(echoRemovedCount) + " echo nodes. " + '\n')
runLog.write ('\n' + "Removed " + str(echoRemovedCount) + " echo nodes. " + '\n' + '\n')
echoList.close()

# Recalculate basic Echo graph statistics
print ("Recalculating various things for Echo graph..." + '\n')
echoNodes = nx.number_of_nodes(echoNW)
echoEdges = nx.number_of_edges(echoNW)
echoDensity = nx.density(echoNW)
echoNodeList = nx.nodes(echoNW)
echoNodeList.sort()
echoSelfLoopTotal = echoNW.number_of_selfloops()
echoConnections = echoEdges - echoSelfLoopTotal

print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Self-loops: ' + str(echoSelfLoopTotal))
print ('Echo Connections (edges minus self-loops): ' + str(echoConnections))
print ('Echo Density: ' + str(echoDensity))
print
print (str(echoNodeList))
print

anFile.write ('\n' + 'New Network Properties: ' + '\n' + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Self-loops: ' + str(echoSelfLoopTotal) + '\n')
anFile.write ('Echo Connections (edges minus self-loops): ' + str(echoConnections) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n' + '\n')

runLog.write ('Echo nodes: ' + '\n' + '\n')
runLog.write (str(echoNodeList) + '\n' + '\n')

# Remove echo self-loops
echoSelfLoopCount = 0
if selfLoopIP == 1:
	print ('Checking for and removing echo self-loops...' + '\n')
	runLog.write('\n' + 'Checking for and removing echo self-loops...' + '\n')

	for u,v in echoNW.edges():
		if u == v:
			echoNW.remove_edge(u,v)
			print ('removed self-loop ' + str(u))
			echoSelfLoopCount += 1

	if echoSelfLoopCount == 0:
		print ('No echo self-loops removed.' + '\n')
		runLog.write('\n' + 'No echo self-loops removed.' + '\n')
	
	if echoSelfLoopCount == 1:
		print ('\n' + 'Removed ' + str(echoSelfLoopCount) + ' echo self-loop edge.')
		runLog.write ('\n' + 'Removed ' + str(echoSelfLoopCount) + ' echo self-loop edge.' + '\n')

	if echoSelfLoopCount > 1:
		print ('\n' + 'Removed ' + str(echoSelfLoopCount) + ' echo self-loop edges.')
		runLog.write ('\n' + 'Removed ' + str(echoSelfLoopCount) + ' echo self-loop edges.' + '\n')

else:
	print ('No echo self-loops removed.' + '\n')
	runLog.write('\n' + 'No echo self-loops removed.' + '\n')

# Count echo isolates
echoIsolateCount = 0
for i in echoNodeList:
	if nx.is_isolate(echoNW,i):
		runLog.write('\n' + 'Echo node ' + i + ' is an Isolate. ' + '\n')
		echoIsolateCount += 1

# Count echo sources and sinks
echoSourceCount = 0
echoSinkCount = 0
if nx.is_directed(echoNW):
	for i in echoNodeList:
		echoOutDeg = echoNW.out_degree(i)
		echoInDeg = echoNW.in_degree(i)
		if echoOutDeg == 0 and echoInDeg >= 1: 
			echoSourceCount += 1

		elif echoInDeg == 0 and echoOutDeg >= 1:
			echoSinkCount += 1
else:
	print("Cannot count sources and sinks. Echo graph is undirected. ")
	print

# Read the MB GEXF and generate graph
print
print ("Reading MB GEXF file and generating graph... ")
mbInputPath = os.path.join("gexf", omegaYear + '_mb.gexf')
mbNW = nx.read_gexf(mbInputPath)

# Calculate basic MB graph statistics
print
print ('Calculating various things for MB graph...' + '\n')
mbNodes = nx.number_of_nodes(mbNW)
mbEdges = nx.number_of_edges(mbNW)
mbDensity = nx.density(mbNW)
mbNodeList = nx.nodes(mbNW)
mbNodeList.sort()
mbSelfLoopTotal = mbNW.number_of_selfloops()
mbConnections = mbEdges - mbSelfLoopTotal

print ("Omega Year: " + omegaYear)
print ('MB Nodes: ' + str(mbNodes))
print ('MB Edges: ' + str(mbEdges))
print ('MB Self-loops: ' + str(mbSelfLoopTotal))
print ('MB Connections (edges minus self-loops): ' + str(mbConnections))
print ('MB Density: ' + str(mbDensity))
print
print (str(nx.info(mbNW)))
print
print (str(mbNodeList))

anFile.write ('Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write ('MB Nodes: ' + str(mbNodes) + '\n')
anFile.write ('MB Edges: ' + str(mbEdges) + '\n')
anFile.write ('MB Self-loops: ' + str(mbSelfLoopTotal) + '\n')
anFile.write ('MB Connections (edges minus self-loops): ' + str(mbConnections) + '\n')
anFile.write ('MB Density: ' + str(mbDensity) + '\n' + '\n')
anFile.write (str(nx.info(mbNW)) + '\n' + '\n')

runLog.write ('\n' + 'MB nodes: ' + '\n' + '\n')
runLog.write (str(mbNodeList) + '\n' + '\n')

# Remove nodes that are not on inclusion-list (`inclusion-lists/mb_inclusion_list.txt')

# Open mb_inclusion_list
mbRemovedCount = 0
print ('\n' + 'Checking list and removing MB nodes...' + '\n')
runLog.write ('Checking list and removing MB nodes...' + '\n' + '\n')
anFile.write ('Checking list and removing MB nodes...' + '\n' + '\n')
	
mbListPath = os.path.join("inclusion-lists", 'master_mb_inclusion_list.txt')
mbList = open(mbListPath, 'r')

# Check for genres and remove if not present
mbIncludeNodes = [omegaYear + '_' + line.strip() for line in mbList]

for i in mbNodeList:
	if not i in mbIncludeNodes:
		mbNW.remove_node(i) 
		print ('Removed MB node ' + str(i))
		runLog.write ('Removed MB node ' + str(i) + '\n')
		# anFile.write (str(i) + '\n')
		mbRemovedCount += 1
	else:
		mbNodeOP = str(i).replace(omegaYear + "_","")
		mbListFile.write(mbNodeOP + '\n')
		#wikiNodeOP = mbwdMapDict[mbNodeOP]
		#wikiListFile.write(wikiNodeOP + '\n')

#wikiListFile.close()

print ('\n'+ "Removed " + str(mbRemovedCount) + " MB nodes. " + '\n')
runLog.write ('\n' + "Removed " + str(mbRemovedCount) + " MB nodes. " + '\n' + '\n')
mbList.close()

# Recalculate basic MB graph statistics
print ("Recalculating various things for MB graph..." + '\n')
mbNodes = nx.number_of_nodes(mbNW)
mbEdges = nx.number_of_edges(mbNW)
mbDensity = nx.density(mbNW)
mbNodeList = nx.nodes(mbNW)
mbNodeList.sort()
mbSelfLoopTotal = mbNW.number_of_selfloops()
mbConnections = mbEdges - mbSelfLoopTotal

print ('MB Nodes: ' + str(mbNodes))
print ('MB Edges: ' + str(mbEdges))
print ('MB Self-loops: ' + str(mbSelfLoopTotal))
print ('MB Connections (edges minus self-loops): ' + str(mbConnections))
print ('MB Density: ' + str(mbDensity))
print
print (str(mbNodeList))
print

anFile.write ('\n' + 'New Network Properties: ' + '\n' + '\n')
anFile.write ('MB Nodes: ' + str(mbNodes) + '\n')
anFile.write ('MB Edges: ' + str(mbEdges) + '\n')
anFile.write ('MB Self-loops: ' + str(mbSelfLoopTotal) + '\n')
anFile.write ('MB Connections (edges minus self-loops): ' + str(mbConnections) + '\n')
anFile.write ('MB Density: ' + str(mbDensity) + '\n' + '\n')

runLog.write ('MB nodes: ' + '\n' + '\n')
runLog.write (str(mbNodeList) + '\n' + '\n')

# Remove MB self-loops
mbSelfLoopCount = 0
if selfLoopIP == 1:
	print ('Checking for and removing MB self-loops...' + '\n')
	runLog.write('\n' + 'Checking for and removing MB self-loops...' + '\n')

	for u,v in mbNW.edges():
		if u == v:
			mbNW.remove_edge(u,v)
			print ('removed self-loop ' + str(u))
			mbSelfLoopCount += 1

	if mbSelfLoopCount == 0:
		print ('No MB self-loops removed.' + '\n')
		runLog.write('\n' + 'No MB self-loops removed.' + '\n')
	
	if mbSelfLoopCount == 1:
		print ('\n' + 'Removed ' + str(mbSelfLoopCount) + ' MB self-loop edge.')
		runLog.write ('\n' + 'Removed ' + str(mbSelfLoopCount) + ' MB self-loop edge.' + '\n')

	if mbSelfLoopCount > 1:
		print ('\n' + 'Removed ' + str(mbSelfLoopCount) + ' MB self-loop edges.')
		runLog.write ('\n' + 'Removed ' + str(mbSelfLoopCount) + ' MB self-loop edges.' + '\n')

else:
	print ('No MB self-loops removed.' + '\n')
	runLog.write('\n' + 'No MB self-loops removed.' + '\n')

# Count MB isolates
mbIsolateCount = 0
for i in mbNodeList:
	if nx.is_isolate(mbNW,i):
		runLog.write('\n' + 'MB node ' + i + ' is an Isolate. ' + '\n')
		mbIsolateCount += 1

# Count MB sources and sinks
mbSourceCount = 0
mbSinkCount = 0
if nx.is_directed(mbNW):
	for i in mbNodeList:
		mbOutDeg = mbNW.out_degree(i)
		mbInDeg = mbNW.in_degree(i)
		if mbOutDeg == 0 and mbInDeg >= 1: 
			mbSourceCount += 1

		elif mbInDeg == 0 and mbOutDeg >= 1:
			mbSinkCount += 1
else:
	print("Cannot count sources and sinks. MB graph is undirected. ")
	print

# Read the Wikidata GEXF and generate graph
print ('\n' + "Reading wikidata GEXF file and generating graph... ")
wdInputPath = os.path.join("gexf", 'wikidata.gexf')
wikiNW = nx.read_gexf(wdInputPath)

# Calculate basic Wiki graph statistics
print ('\n' + 'Calculating various things for Wiki graph...' + '\n')
wikiNodes = nx.number_of_nodes(wikiNW)
wikiEdges = nx.number_of_edges(wikiNW)
wikiDensity = nx.density(wikiNW)
wikiNodeList = nx.nodes(wikiNW)
wikiNodeList.sort()
wikiSelfLoopTotal = wikiNW.number_of_selfloops()
wikiConnections = wikiEdges - wikiSelfLoopTotal

print ("Omega Year: " + omegaYear)
print ('Wiki Nodes: ' + str(wikiNodes))
print ('Wiki Edges: ' + str(wikiEdges))
print ('Wiki Self-loops: ' + str(wikiSelfLoopTotal))
print ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections))
print ('Wiki Density: ' + str(wikiDensity))
print
print (str(nx.info(wikiNW)))
print
print (str(wikiNodeList))

anFile.write ('Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')
anFile.write (str(nx.info(wikiNW)) + '\n' + '\n')

runLog.write ("Omega Year: " + omegaYear + '\n' + '\n')
runLog.write ('Wiki nodes: ' + '\n' + '\n')
runLog.write (str(wikiNodeList) + '\n' + '\n')

# Remove nodes that are not on inclusion-lists (`inclusion-lists/wiki_inclusion_list.txt' and `inclusion-lists/echo_inclusion_list.txt')
# Open wiki_inclusion_list

wikiRemovedCount = 0
print ('\n' + 'Checking list and removing wiki nodes...' + '\n')
runLog.write ('Checking list and removing wiki nodes...' + '\n' + '\n')
anFile.write ('Checking list and removing wiki nodes...' + '\n' + '\n')

wikiListPath = os.path.join("inclusion-lists", omegaYear + '_wiki_inclusion_list.txt')
wikiList = open(wikiListPath, 'r')
wikiIncludeNodes = [line.strip() for line in wikiList]

for i in wikiNodeList:
	if not i in wikiIncludeNodes:
		wikiNW.remove_node(i) 
		print ('Removed wiki node ' + str(i))
		runLog.write ('Removed wiki node ' + str(i) + '\n')
		# anFile.write (str(i) + '\n')
		wikiRemovedCount += 1

print ('\n'+ "Removed " + str(wikiRemovedCount) + " wiki nodes. " + '\n')
runLog.write ('\n' + "Removed " + str(wikiRemovedCount) + " wiki nodes. " + '\n' + '\n')
wikiList.close()

# Recalculate basic Wiki graph statistics
print ('Recalculating various things for Wiki graph...' + '\n')
wikiNodes = nx.number_of_nodes(wikiNW)
wikiEdges = nx.number_of_edges(wikiNW)
wikiDensity = nx.density(wikiNW)
wikiNodeList = nx.nodes(wikiNW)
wikiNodeList.sort()
wikiSelfLoopTotal = wikiNW.number_of_selfloops()
wikiConnections = wikiEdges - wikiSelfLoopTotal

print ('Wiki Nodes: ' + str(wikiNodes))
print ('Wiki Edges: ' + str(wikiEdges))
print ('Wiki Self-loops: ' + str(wikiSelfLoopTotal))
print ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections))
print ('Wiki Density: ' + str(wikiDensity))
print
print (str(wikiNodeList))
print

anFile.write ('\n' + '\n' + 'New Network Properties: ' + '\n' + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')

runLog.write ('Wiki nodes: ' + '\n' + '\n')
runLog.write (str(wikiNodeList) + '\n' + '\n')

# Remove wiki self-loops
wikiSelfLoopCount = 0
if selfLoopIP == 1:
	print ('Checking for and removing wiki self-loops...' + '\n')
	runLog.write('Checking for and removing wiki self-loops...' + '\n')

	for u,v in wikiNW.edges():
		if u == v:
			wikiNW.remove_edge(u,v)
			print ('removed self-loop ' + str(u))
			wikiSelfLoopCount += 1

	if wikiSelfLoopCount == 0:
		print ('No wiki self-loops removed.' + '\n')
		runLog.write('\n' + 'No wiki self-loops removed.' + '\n')
	
	if wikiSelfLoopCount == 1:
		print ('\n' + 'Removed ' + str(wikiSelfLoopCount) + ' wiki self-loop edge.')
		runLog.write ('\n' + 'Removed ' + str(wikiSelfLoopCount) + ' wiki self-loop edge.' + '\n')

	if wikiSelfLoopCount > 1:
		print ('\n' + 'Removed ' + str(wikiSelfLoopCount) + ' wiki self-loop edges.')
		runLog.write ('\n' + 'Removed ' + str(wikiSelfLoopCount) + ' wiki self-loop edges.' + '\n')

else:
	print ('No wiki self-loops removed.' + '\n')
	runLog.write('\n' + 'No wiki self-loops removed.' + '\n')

# Count wiki isolates
wikiIsolateCount = 0
for i in wikiNodeList:
	if nx.is_isolate(wikiNW,i):
		runLog.write('\n' + 'Wiki node ' + i + ' is an Isolate. ' + '\n')
		wikiIsolateCount += 1

# Count wiki sources and sinks
wikiSourceCount = 0
wikiSinkCount = 0
if nx.is_directed(wikiNW):
	for i in wikiNodeList:
		wikiOutDeg = wikiNW.out_degree(i)
		wikiInDeg = wikiNW.in_degree(i)
		if wikiOutDeg == 0 and wikiInDeg >= 1: 
			wikiSourceCount += 1

		elif wikiInDeg == 0 and wikiOutDeg >= 1:
			wikiSinkCount += 1
else:
	print("Cannot count sources and sinks. Wiki graph is undirected. ")
	print

# Recalculate basic Wiki graph statistics
print ('Recalculating various things for Wiki graph...' + '\n')
wikiNodes = nx.number_of_nodes(wikiNW)
wikiEdges = nx.number_of_edges(wikiNW)
wikiDensity = nx.density(wikiNW)
wikiNodeList = nx.nodes(wikiNW)
wikiNodeList.sort()
wikiSelfLoopTotal = wikiNW.number_of_selfloops()
wikiConnections = wikiEdges - wikiSelfLoopTotal

print ('Wiki Nodes: ' + str(wikiNodes))
print ('Wiki Isolates: ' + str(wikiIsolateCount))
print ('Wiki Edges: ' + str(wikiEdges))
print ('Wiki Self-loops: ' + str(wikiSelfLoopTotal))
print ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections))
print ('Wiki Density: ' + str(wikiDensity))
print

anFile.write ('Final Network Properties: ' + '\n' + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Isolates: ' + str(wikiIsolateCount) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')

# Write gexf files for use in Gephi
print
print ("Writing gexf files for use in Gephi... " + '\n')
runLog.write('\n' + "Writing gexf files for use in Gephi... " + '\n')
nx.write_gexf(wikiNW, wikiGexfFile)
wikiGexfFile.close()
nx.write_gexf(echoNW, echoGexfFile)
echoGexfFile.close()
nx.write_gexf(mbNW, mbGexfFile)
mbGexfFile.close()

# Plot and display graphs
# Graph plotting parameters - moved to config file 'config_nw.txt'
print ('Reading layout config file... ' + '\n')
runLog.write('\n' + "Reading layout config file... " + '\n')

# Open and read 'config_nw.txt'
nwConfigPath = os.path.join ("config", 'config_nw.txt')
nwConfig = open(nwConfigPath, 'r').readlines()

# Remove the first line
firstLine = nwConfig.pop(0)

for line in nwConfig:
	n_size, n_alpha, node_colour, n_text_size, text_font, e_thickness, e_alpha, edge_colour, l_pos, e_text_size, edge_label_colour = line.split(",")
	
node_size = int(n_size)
node_alpha = float(n_alpha)
node_text_size = int(n_text_size)
edge_thickness = int(e_thickness)
edge_alpha = float(e_alpha)
label_pos = float(l_pos)
edge_text_size = int(e_text_size)

print ('Laying out wiki graph... ' + '\n')
runLog.write('\n' + "Laying out wiki graph... " + '\n')

graph_pos = nx.spring_layout(wikiNW)
nx.draw_networkx_nodes(wikiNW, graph_pos, node_size = node_size, alpha = node_alpha, node_color=node_colour)
nx.draw_networkx_edges(wikiNW, graph_pos, width = edge_thickness, alpha = edge_alpha, color = edge_colour)
#nx.draw_networkx_labels(wikiNW, graph_pos, font_size = node_text_size, font_family = text_font)

# write image file
print ('Writing wiki network image file... ' + '\n')
runLog.write('\n' + "Writing wiki network image file... " + '\n')

plt.savefig(wikiImgPath, format = 'eps', bbox_inches='tight')
#wikiImg.close()

# display graph
print ('Displaying wiki graph...' + '\n')
plt.show()

# Clear plot
plt.clf()

print ('Laying out echo graph...' + '\n')
runLog.write('\n' + "Laying out echo graph... " + '\n')

graph_pos = nx.spring_layout(echoNW)
nx.draw_networkx_nodes(echoNW, graph_pos, node_size = node_size, alpha = node_alpha, node_color=node_colour)
nx.draw_networkx_edges(echoNW, graph_pos, width = edge_thickness, alpha = edge_alpha, color = edge_colour)
#nx.draw_networkx_labels(echoNW, graph_pos, font_size = node_text_size, font_family = text_font)

# write image file
print ('Writing echo network image file...' + '\n')
runLog.write('\n' + "Writing echo network image file... " + '\n')

plt.savefig(echoImgPath, format = 'eps', bbox_inches='tight')
#echoImg.close()

# display graph
print ('Displaying echo graph...' + '\n')
plt.show()

# Clear plot
plt.clf()

print ('Laying out MB graph...' + '\n')
runLog.write('\n' + "Laying out MB graph... " + '\n')

graph_pos = nx.spring_layout(mbNW)
nx.draw_networkx_nodes(mbNW, graph_pos, node_size = node_size, alpha = node_alpha, node_color=node_colour)
nx.draw_networkx_edges(mbNW, graph_pos, width = edge_thickness, alpha = edge_alpha, color = edge_colour)
#nx.draw_networkx_labels(echoNW, graph_pos, font_size = node_text_size, font_family = text_font)

# write image file
print ('Writing MB network image file...' + '\n')
runLog.write('\n' + "Writing MB network image file... " + '\n')

plt.savefig(mbImgPath, format = 'eps', bbox_inches='tight')
#mbImg.close()

# display graph
print ('Displaying MB graph...' + '\n')
plt.show()

# Clear plot
plt.clf()

# Analysis of Wiki graph
print
print ('Analysing Wiki graph...' + '\n')
wikiNodes = nx.number_of_nodes(wikiNW)
wikiEdges = nx.number_of_edges(wikiNW)
wikiDensity = nx.density(wikiNW)
wikiDegreeSeq = sorted(nx.degree(wikiNW).values(),reverse=True)
wikiMaxDeg = max(wikiDegreeSeq)

try:
	if nx.is_directed(wikiNW):
		wikiAvClustering = 0
		wikiConnectComp = 0
		wikiIsDigraph = ("Directed")
		wikiUndNW = wikiNW.to_undirected()
		wikiPart = community.best_partition(wikiUndNW)
		wikiMod = community.modularity(wikiPart, wikiUndNW)
		wikiAvClustering = nx.average_clustering(wikiUndNW)
		wikiConnectComp = [len(c) for c in sorted(nx.connected_components(wikiUndNW), key=len, reverse=True)]
	else:
		wikiPart = community.best_partition(wikiNW)
		wikiMod = community.modularity(wikiPart, wikiNW)
		wikiAvClustering = nx.average_clustering(wikiNW)
		wikiConnectComp = [len(c) for c in sorted(nx.connected_components(wikiNW), key=len, reverse=True)]
		wikiIsDigraph = ("Undirected")
except:
	wikiMod = 0
	wikiAvClustering = 0
	wikiConnectComp = 0

print ("Wiki graph is " + wikiIsDigraph)
print ("Wiki Nodes: " + str(wikiNodes))
print ("Wiki Isolates: " + str(wikiIsolateCount))
print ("Wiki Edges: " + str(wikiEdges))
print ("Wiki Density: " + str(wikiDensity))
print ("Wiki Maximal degree: " + str(wikiMaxDeg))
print ("Wiki modularity: " + str(wikiMod))
print ('Wiki Average Clustering Coefficient: ' + str(wikiAvClustering))
print ('Wiki Connected Components: ' + str(wikiConnectComp))
print ("Source nodes: " + str(wikiSourceCount))
print ("Sink nodes: " + str(wikiSinkCount))
print

anFile.write ('\n' + 'Wiki Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write ("Wiki graph is " + wikiIsDigraph + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Isolates: ' + str(wikiIsolateCount) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n')
anFile.write ('Wiki Maximal degree: ' + str(wikiMaxDeg) + '\n')
anFile.write ("Wiki modularity: " + str(wikiMod) + '\n')
anFile.write ('Wiki Average Clustering Coefficient: ' + str(wikiAvClustering) + '\n')
anFile.write ('Wiki Connected Components: ' + str(wikiConnectComp) + '\n')
anFile.write ("Source nodes: " + str(wikiSourceCount) + '\n')
anFile.write ("Sink nodes: " + str(wikiSinkCount) + '\n' + '\n')

try:
	for key,value in wikiPart.items():
		wikiCommFile.write(key + "," + str(value) + '\n')
except:
	pass

anFile.write ('\n' + 'Confirmation: ' + '\n' + str(nx.info(wikiNW)) + '\n')

runLog.write ("Omega Year: " + omegaYear + '\n' + '\n')
runLog.write('\n' + "Wiki graph is " + wikiIsDigraph + '\n')
runLog.write ('\n' + "Properties: " + '\n' + str(nx.info(wikiNW)) + '\n')

# Recalculate basic Echo graph statistics
print ('Analysing Echo graph...' + '\n')
echoNodes = nx.number_of_nodes(echoNW)
echoEdges = nx.number_of_edges(echoNW)
echoDensity = nx.density(echoNW)
echoNodeList = nx.nodes(echoNW)
echoNodeList.sort()
echoSelfLoopTotal = echoNW.number_of_selfloops()
echoConnections = echoEdges - echoSelfLoopTotal
echoDegreeSeq = sorted(nx.degree(echoNW).values(),reverse=True)
echoMaxDeg = max(echoDegreeSeq)

print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Self-loops: ' + str(echoSelfLoopTotal))
print ('Echo Connections (edges minus self-loops): ' + str(echoConnections))
print ('Echo Density: ' + str(echoDensity))
print ('Echo Isolates: ' + str(echoIsolateCount))
print

anFile.write ('\n' + 'Final Network Properties: ' + '\n' + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Self-loops: ' + str(echoSelfLoopTotal) + '\n')
anFile.write ('Echo Connections (edges minus self-loops): ' + str(echoConnections) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n')
anFile.write ('Echo Isolates: ' + str(echoIsolateCount) + '\n' + '\n')

# Calculate modularity/community memberships
try:
	if nx.is_directed(echoNW):
		echoAvClustering = 0
		echoConnectComp = 0
		echoIsDigraph = ("Directed")
		echoUndNW = echoNW.to_undirected()
		echoPart = community.best_partition(echoUndNW)
		echoMod = community.modularity(echoPart, echoUndNW)
		echoAvClustering = nx.average_clustering(echoUndNW)
		echoConnectComp = [len(c) for c in sorted(nx.connected_components(echoUndNW), key=len, reverse=True)]
	else:
		echoPart = community.best_partition(echoNW)
		echoMod = community.modularity(echoPart, echoNW)
		echoAvClustering = nx.average_clustering(echoNW)
		echoConnectComp = [len(c) for c in sorted(nx.connected_components(echoNW), key=len, reverse=True)]
		echoIsDigraph = ("Undirected")
except:
	echoMod = 0
	echoAvClustering = 0
	echoConnectComp = 0

print ("Echo graph is " + echoIsDigraph)
print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Isolates: ' + str(echoIsolateCount))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Density: ' + str(echoDensity))
print ('Echo Maximal degree: ' + str(echoMaxDeg))
print ("Echo modularity: " + str(echoMod))
print ('Echo Average Clustering Coefficient: ' + str(echoAvClustering))
print ('Echo Connected Components: ' + str(echoConnectComp))
print ("Source nodes: " + str(echoSourceCount))
print ("Sink nodes: " + str(echoSinkCount))
print

anFile.write ('\n' + 'Echo Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write("Echo graph is " + echoIsDigraph + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Isolates: ' + str(echoIsolateCount) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n')
anFile.write ('Echo Maximal degree: ' + str(echoMaxDeg) + '\n')
anFile.write ("Echo modularity: " + str(echoMod) + '\n')
anFile.write ('Echo Average Clustering Coefficient: ' + str(echoAvClustering) + '\n')
anFile.write ('Echo Connected Components: ' + str(echoConnectComp) + '\n')
anFile.write ("Source nodes: " + str(echoSourceCount) + '\n')
anFile.write ("Sink nodes: " + str(echoSinkCount) + '\n' + '\n')

try:
	for key,value in echoPart.items():
		genLabel = key.replace("2015_", omegaYear + '_')
		echoCommFile.write(genLabel + "," + str(value) + '\n')
except:
	pass

anFile.write ('\n' + 'Confirmation: ' + '\n' + str(nx.info(echoNW)))

runLog.write('\n' + "Echo graph is " + echoIsDigraph + '\n')
runLog.write ('\n' + "Properties: " + '\n' + str(nx.info(echoNW)) + '\n')

# Recalculate basic MB graph statistics
print ('Analysing MB graph...' + '\n')
mbNodes = nx.number_of_nodes(mbNW)
mbEdges = nx.number_of_edges(mbNW)
mbDensity = nx.density(mbNW)
mbNodeList = nx.nodes(mbNW)
mbNodeList.sort()
mbSelfLoopTotal = mbNW.number_of_selfloops()
mbConnections = mbEdges - mbSelfLoopTotal
mbDegreeSeq = sorted(nx.degree(mbNW).values(),reverse=True)
mbMaxDeg = max(mbDegreeSeq)

print ('MB Nodes: ' + str(mbNodes))
print ('MB Edges: ' + str(mbEdges))
print ('MB Self-loops: ' + str(mbSelfLoopTotal))
print ('MB Connections (edges minus self-loops): ' + str(mbConnections))
print ('MB Density: ' + str(mbDensity))
print ('MB Isolates: ' + str(mbIsolateCount))
print

anFile.write ('\n' + 'Final Network Properties: ' + '\n' + '\n')
anFile.write ('MB Nodes: ' + str(mbNodes) + '\n')
anFile.write ('MB Edges: ' + str(mbEdges) + '\n')
anFile.write ('MB Self-loops: ' + str(mbSelfLoopTotal) + '\n')
anFile.write ('MB Connections (edges minus self-loops): ' + str(mbConnections) + '\n')
anFile.write ('MB Density: ' + str(mbDensity) + '\n')
anFile.write ('MB Isolates: ' + str(mbIsolateCount) + '\n' + '\n')

# Calculate modularity/community memberships
try:
	if nx.is_directed(mbNW):
		mbAvClustering = 0
		mbConnectComp = 0
		mbIsDigraph = ("Directed")
		mbUndNW = mbNW.to_undirected()
		mbPart = community.best_partition(mbUndNW)
		mbMod = community.modularity(mbPart, mbUndNW)
		mbAvClustering = nx.average_clustering(mbUndNW)
		mbConnectComp = [len(c) for c in sorted(nx.connected_components(mbUndNW), key=len, reverse=True)]
	else:
		mbPart = community.best_partition(mbNW)
		mbMod = community.modularity(mbPart, mbNW)
		mbAvClustering = nx.average_clustering(mbNW)
		mbConnectComp = [len(c) for c in sorted(nx.connected_components(mbNW), key=len, reverse=True)]
		mbIsDigraph = ("Undirected")
except:
	mbMod = 0
	mbAvClustering = 0
	mbConnectComp = 0

print ("MB graph is " + mbIsDigraph)
print ('MB Nodes: ' + str(mbNodes))
print ('MB Isolates: ' + str(mbIsolateCount))
print ('MB Edges: ' + str(mbEdges))
print ('MB Density: ' + str(mbDensity))
print ('MB Maximal degree: ' + str(mbMaxDeg))
print ("MB modularity: " + str(mbMod))
print ('MB Average Clustering Coefficient: ' + str(mbAvClustering))
print ('MB Connected Components: ' + str(mbConnectComp))
print ("Source nodes: " + str(mbSourceCount))
print ("Sink nodes: " + str(mbSinkCount))
print

anFile.write ('\n' + 'MB Network Properties: ' + '\n' + '\n')
anFile.write ("Omega Year: " + omegaYear + '\n')
anFile.write("MB graph is " + mbIsDigraph + '\n')
anFile.write ('MB Nodes: ' + str(mbNodes) + '\n')
anFile.write ('MB Isolates: ' + str(mbIsolateCount) + '\n')
anFile.write ('MB Edges: ' + str(mbEdges) + '\n')
anFile.write ('MB Density: ' + str(mbDensity) + '\n')
anFile.write ('MB Maximal degree: ' + str(mbMaxDeg) + '\n')
anFile.write ("MB modularity: " + str(mbMod) + '\n')
anFile.write ('MB Average Clustering Coefficient: ' + str(mbAvClustering) + '\n')
anFile.write ('MB Connected Components: ' + str(mbConnectComp) + '\n')
anFile.write ("Source nodes: " + str(mbSourceCount) + '\n')
anFile.write ("Sink nodes: " + str(mbSinkCount) + '\n' + '\n')

try:
	for key,value in mbPart.items():
		genLabel = key.replace("2015_", omegaYear + '_')
		mbCommFile.write(genLabel + "," + str(value) + '\n')
except:
	pass

anFile.write ('\n' + 'Confirmation: ' + '\n' + str(nx.info(mbNW)))
anFile.close()

runLog.write('\n' + "MB graph is " + mbIsDigraph + '\n')
runLog.write ('\n' + "Properties: " + '\n' + str(nx.info(mbNW)) + '\n')

# End timing of run
endTime = datetime.now()

print
print ("Omega Year: " + omegaYear + '\n')
print ('Date of run: {}'.format(runDate))
print ('Duration of run : {}'.format(endTime - startTime))

runLog.write ("Omega Year: " + omegaYear + '\n' + '\n')
runLog.write ('\n' + '\n' + 'Date of run: {}'.format(runDate) + '\n')
runLog.write ('Duration of run : {}'.format(endTime - startTime) + '\n')
runLog.close()
