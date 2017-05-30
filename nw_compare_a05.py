# nw_compare_a05.py
# Version a05
# by jmg - j.gagen*AT*gold*DOT*ac*DOT*uk
# May 30th 2017

# Licence: http://creativecommons.org/licenses/by-nc-sa/3.0/

# Builds networks from input GEXFs, then removes nodes based on inclusion-lists
# Builds both EchoNest and Wikidata networks and compares
# This version applies different analyses for undirected and directed graphs
# Writes modularity/community files to `results/communities' when graph is undirected

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

versionNumber = ("a05")

# Initiate timing of run
runDate = datetime.now()
startTime = datetime.now()

# Create 'logs' subdirectory if necessary
if not os.path.exists("logs"):
    os.makedirs("logs")

# Create 'gexf' subdirectories if necessary
if not os.path.exists("gexf"):
	os.makedirs("gexf")

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

# Open file for writing log
logPath = os.path.join("logs", 'nw_compare_' + versionNumber  + '_' + str(runDate) + '_' + str(startTime) + '_log.txt')
runLog = open(logPath, 'a')

# Begin
print ('\n' + "Network Comparison Thing | Version " + versionNumber + " | Starting...")
runLog.write ("==========================================================================" + '\n' + '\n')
runLog.write ("Network Comparison Thing | Version " + versionNumber + '\n' + '\n')

# Open files for writing gexfs
wikiGexfPath = os.path.join("gexf", 'nw_compare_wiki_' + versionNumber + '.gexf')
wikiGexfFile = open(wikiGexfPath, 'w')

echoGexfPath = os.path.join("gexf", 'nw_compare_echo_' + versionNumber + '.gexf')
echoGexfFile = open(echoGexfPath, 'w')

# Open file for `communities' ouput files
wikiCommPath = os.path.join("results/communities", 'nw_compare_wiki_' + versionNumber + '_' + str(runDate) + '_' + str(startTime) + '_communities.txt')
wikiCommFile = open(wikiCommPath, 'w')

echoCommPath = os.path.join("results/communities", 'nw_compare_echo_' + versionNumber + '_' + str(runDate) + '_' + str(startTime) + '_communities.txt')
echoCommFile = open(echoCommPath, 'w')

# Open file for analysis results
anPath = os.path.join("results/analysis", 'nw_compare_' + versionNumber  + '_' + str(runDate) + '_' + str(startTime) + '_analysis.txt')
anFile = open(anPath, 'w')

# Open files to write images
wikiImgPath = os.path.join("networks", 'nw_compare_wiki_' + versionNumber + '_nw.eps')
wikiImg = open (wikiImgPath, 'w')

echoImgPath = os.path.join("networks", 'nw_compare_echo_' + versionNumber + '_nw.eps')
echoImg = open (echoImgPath, 'w')

anFile.write ("==========================================================================" + '\n' + '\n')
anFile.write ("Network Comparison Thing | Version " + versionNumber + '\n' + '\n')

# Uncomment the line below to facilitate optional self-loop removal
#selfLoopIP = int(input ("Enter 1 here to remove self-loop edges: "))
# Auto-remove self-loops - comment out the line below if self-loop removal is optional
selfLoopIP = 1

# Read the Wikidata GEXF and generate graph
print ('\n' + "Reading wikidata GEXF file and generating graph... ")
wdInputPath = os.path.join("gexf", 'wikidata.gexf')
wikiNW = nx.read_gexf(wdInputPath)

# Read the EchiNest GEXF and generate graph
print ('\n' + "Reading echoNest GEXF file and generating graph... ")
enInputPath = os.path.join("gexf", 'echonest.gexf')
echoNW = nx.read_gexf(enInputPath)

# Calculate basic Wiki graph statistics
print ('\n' + 'Calculating various things for Wiki graph...' + '\n')
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
print (str(nx.info(wikiNW)))
print
print (str(wikiNodeList))

anFile.write ('Network Properties: ' + '\n' + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')
anFile.write (str(nx.info(wikiNW)) + '\n' + '\n')

runLog.write ('Wiki nodes: ' + '\n' + '\n')
runLog.write (str(wikiNodeList) + '\n' + '\n')

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
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Self-loops: ' + str(echoSelfLoopTotal) + '\n')
anFile.write ('Echo Connections (edges minus self-loops): ' + str(echoConnections) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n' + '\n')
anFile.write (str(nx.info(echoNW)) + '\n' + '\n')

runLog.write ('Echo nodes: ' + '\n' + '\n')
runLog.write (str(echoNodeList) + '\n' + '\n')

# Remove nodes that are not on inclusion-lists (`inclusion-lists/wiki_inclusion_list.txt' and `inclusion-lists/echo_inclusion_list.txt')
# Open wiki_inclusion_list
wikiRemovedCount = 0
print ('\n' + 'Checking list and removing wiki nodes...' + '\n')
runLog.write ('Checking list and removing wiki nodes...' + '\n' + '\n')
anFile.write ('Checking list and removing wiki nodes...' + '\n' + '\n')

wikiListPath = os.path.join("inclusion-lists", 'wiki_inclusion_list.txt')
wikiList = open(wikiListPath, 'r')
wikiIncludeNodes = [line.strip() for line in wikiList]

for i in wikiNodeList:
	if not i in wikiIncludeNodes:
		wikiNW.remove_node(i) 
		print ('Removed wiki node ' + str(i))
		runLog.write ('Removed wiki node ' + str(i) + '\n')
		wikiRemovedCount += 1

print ('\n'+ "Removed " + str(wikiRemovedCount) + " wiki nodes. " + '\n')
runLog.write ('\n' + "Removed " + str(wikiRemovedCount) + " wiki nodes. " + '\n' + '\n')
wikiList.close()

# Open echo_inclusion_list
echoRemovedCount = 0
print ('\n' + 'Checking list and removing echo nodes...' + '\n')
runLog.write ('Checking list and removing echo nodes...' + '\n' + '\n')
anFile.write ('Checking list and removing echo nodes...' + '\n' + '\n')

echoListPath = os.path.join("inclusion-lists", 'echo_inclusion_list.txt')
echoList = open(echoListPath, 'r')
echoIncludeNodes = [line.strip() for line in echoList]

for i in echoNodeList:
	if not i in echoIncludeNodes:
		echoNW.remove_node(i) 
		print ('Removed echo node ' + str(i))
		runLog.write ('Removed echo node ' + str(i) + '\n')
		echoRemovedCount += 1

print ('\n'+ "Removed " + str(echoRemovedCount) + " echo nodes. " + '\n')
runLog.write ('\n' + "Removed " + str(echoRemovedCount) + " echo nodes. " + '\n' + '\n')
echoList.close()

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

anFile.write ('\n' + 'New Network Properties: ' + '\n' + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')

runLog.write ('Wiki nodes: ' + '\n' + '\n')
runLog.write (str(wikiNodeList) + '\n' + '\n')

# Recalculate basic Echo graph statistics
print ('\n' + "Recalculating various things for Echo graph..." + '\n')
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

anFile.write ('New New Network Properties: ' + '\n' + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Isolates: ' + str(wikiIsolateCount) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Self-loops: ' + str(wikiSelfLoopTotal) + '\n')
anFile.write ('Wiki Connections (edges minus self-loops): ' + str(wikiConnections) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n' + '\n')

# Recalculate basic Echo graph statistics
print ('Recalculating various things for Echo graph...' + '\n')
echoNodes = nx.number_of_nodes(echoNW)
echoEdges = nx.number_of_edges(echoNW)
echoDensity = nx.density(echoNW)
echoNodeList = nx.nodes(echoNW)
echoNodeList.sort()
echoSelfLoopTotal = echoNW.number_of_selfloops()
echoConnections = echoEdges - echoSelfLoopTotal

print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Isolates: ' + str(echoIsolateCount))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Self-loops: ' + str(echoSelfLoopTotal))
print ('Echo Connections (edges minus self-loops): ' + str(echoConnections))
print ('Echo Density: ' + str(echoDensity))
print

anFile.write ('New New Network Properties: ' + '\n' + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Isolates: ' + str(echoIsolateCount) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Self-loops: ' + str(echoSelfLoopTotal) + '\n')
anFile.write ('Echo Connections (edges minus self-loops): ' + str(echoConnections) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n' + '\n')

# Write gexf files for use in Gephi
print
print ("Writing gexf files for use in Gephi... " + '\n')
runLog.write('\n' + "Writing gexf files for use in Gephi... " + '\n')
nx.write_gexf(wikiNW, wikiGexfFile)
wikiGexfFile.close()
nx.write_gexf(echoNW, echoGexfFile)
echoGexfFile.close()

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

plt.savefig(wikiImg, format = 'eps', bbox_inches='tight')
wikiImg.close()

# display graph
print ('Displaying wiki graph...' + '\n')
plt.show()

print ('Laying out echo graph...' + '\n')
runLog.write('\n' + "Laying out echo graph... " + '\n')

graph_pos = nx.spring_layout(echoNW)
nx.draw_networkx_nodes(echoNW, graph_pos, node_size = node_size, alpha = node_alpha, node_color=node_colour)
nx.draw_networkx_edges(echoNW, graph_pos, width = edge_thickness, alpha = edge_alpha, color = edge_colour)
#nx.draw_networkx_labels(echoNW, graph_pos, font_size = node_text_size, font_family = text_font)

# write image file
print ('Writing echo network image file...' + '\n')
runLog.write('\n' + "Writing echo network image file... " + '\n')

plt.savefig(echoImg, format = 'eps', bbox_inches='tight')
echoImg.close()

# display graph
print ('Displaying echo graph...' + '\n')
plt.show()

# Analysis of Wiki graph
print
print ('Analysing Wiki graph...' + '\n')
wikiNodes = nx.number_of_nodes(wikiNW)
wikiEdges = nx.number_of_edges(wikiNW)
wikiDensity = nx.density(wikiNW)
wikiDegreeSeq = sorted(nx.degree(wikiNW).values(),reverse=True)
wikiMaxDeg = max(wikiDegreeSeq)

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

print ("Wiki graph is " + wikiIsDigraph)
print ("Wiki Nodes: " + str(wikiNodes))
print ("Wiki Isolates: " + str(wikiIsolateCount))
print ("Wiki Edges: " + str(wikiEdges))
print ("Wiki Density: " + str(wikiDensity))
print ("Wiki Maximal degree: " + str(wikiMaxDeg))

if wikiAvClustering != 0:
	print ("Wiki modularity: " + str(wikiMod))
	print ('Wiki Average Clustering Coefficient: ' + str(wikiAvClustering))
	print ('Wiki Connected Components: ' + str(wikiConnectComp))
else:
	print ("Source nodes: " + str(wikiSourceCount))
	print ("Sink nodes: " + str(wikiSinkCount))

print

anFile.write ('\n' + 'Wiki Network Properties: ' + '\n' + '\n')
anFile.write("Wiki graph is " + wikiIsDigraph + '\n')
anFile.write ('Wiki Nodes: ' + str(wikiNodes) + '\n')
anFile.write ('Wiki Isolates: ' + str(wikiIsolateCount) + '\n')
anFile.write ('Wiki Edges: ' + str(wikiEdges) + '\n')
anFile.write ('Wiki Density: ' + str(wikiDensity) + '\n')
anFile.write ('Wiki Maximal degree: ' + str(wikiMaxDeg) + '\n')

if wikiAvClustering != 0:
	anFile.write ("Wiki modularity: " + str(wikiMod) + '\n')

	for key,value in wikiPart.items():
		wikiCommFile.write(key + "," + str(value) + '\n')

	anFile.write ('Wiki Average Clustering Coefficient: ' + str(wikiAvClustering) + '\n')
	anFile.write ('Wiki Connected Components: ' + str(wikiConnectComp) + '\n')
else:
	anFile.write ("Source nodes: " + str(wikiSourceCount) + '\n')
	anFile.write ("Sink nodes: " + str(wikiSinkCount) + '\n')

anFile.write ('\n' + 'Confirmation: ' + '\n' + str(nx.info(wikiNW)) + '\n')

runLog.write('\n' + "Wiki graph is " + wikiIsDigraph + '\n')
runLog.write ('\n' + "Properties: " + '\n' + str(nx.info(wikiNW)) + '\n')

# Analysis of Echo graph
print ('Analysing Echo graph...' + '\n')
echoNodes = nx.number_of_nodes(echoNW)
echoEdges = nx.number_of_edges(echoNW)
echoDensity = nx.density(echoNW)
echoDegreeSeq = sorted(nx.degree(echoNW).values(),reverse=True)
echoMaxDeg = max(echoDegreeSeq)

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

print ("Echo graph is " + echoIsDigraph)
print ('Echo Nodes: ' + str(echoNodes))
print ('Echo Isolates: ' + str(echoIsolateCount))
print ('Echo Edges: ' + str(echoEdges))
print ('Echo Density: ' + str(echoDensity))
print ('Echo Maximal degree: ' + str(echoMaxDeg))

if echoAvClustering != 0:
	print ("Echo modularity: " + str(echoMod))
	print ('Echo Average Clustering Coefficient: ' + str(echoAvClustering))
	print ('Echo Connected Components: ' + str(echoConnectComp))
else:
	print ("Source nodes: " + str(echoSourceCount))
	print ("Sink nodes: " + str(echoSinkCount))

print

anFile.write ('\n' + 'Echo Network Properties: ' + '\n' + '\n')
anFile.write("Echo graph is " + echoIsDigraph + '\n')
anFile.write ('Echo Nodes: ' + str(echoNodes) + '\n')
anFile.write ('Echo Isolates: ' + str(echoIsolateCount) + '\n')
anFile.write ('Echo Edges: ' + str(echoEdges) + '\n')
anFile.write ('Echo Density: ' + str(echoDensity) + '\n')
anFile.write ('Echo Maximal degree: ' + str(echoMaxDeg) + '\n')

if echoAvClustering != 0:
	anFile.write ("Echo modularity: " + str(echoMod) + '\n')

	for key,value in echoPart.items():
		echoCommFile.write(key + "," + str(value) + '\n')

	anFile.write ('Echo Average Clustering Coefficient: ' + str(echoAvClustering) + '\n')
	anFile.write ('Echo Connected Components: ' + str(echoConnectComp) + '\n')
else:
	anFile.write ("Source nodes: " + str(echoSourceCount) + '\n')
	anFile.write ("Sink nodes: " + str(echoSinkCount) + '\n')

anFile.write ('\n' + 'Confirmation: ' + '\n' + str(nx.info(echoNW)))
anFile.close()

runLog.write('\n' + "Echo graph is " + echoIsDigraph + '\n')
runLog.write ('\n' + "Properties: " + '\n' + str(nx.info(echoNW)) + '\n')

# End timing of run
endTime = datetime.now()

print
print ('Date of run: {}'.format(runDate))
print ('Duration of run : {}'.format(endTime - startTime))

runLog.write ('\n' + '\n' + 'Date of run: {}'.format(runDate) + '\n')
runLog.write ('Duration of run : {}'.format(endTime - startTime) + '\n')
runLog.close()