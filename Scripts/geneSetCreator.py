# -*- coding: utf-8 -*-
### The working functions

### Joseph Hostyk

### 2018-4-9

from helper import *



###################################################################################################

### First steps: compare our results to GTR's.

###################################################################################################
###
def getKeyTerms():

	GTRfilesPath = PROJECT + "Results/GeneSets/GTRgeneLists/"

	filenames = os.listdir(GTRfilesPath)
	keyTerms = [filename.replace(".tsv", "") for filename in filenames]

	print "We have {} keyterms.".format(len(keyTerms))
	return keyTerms

def createOurGeneSets(keyWords, outputFolder):

	allSynopses = getAllSynopses()
	officialGeneList = getOfficialGenes()
	# inheritanceDictionary = makeDictFromTwoColumns(PROJECT + "allGenes/{}/allGenesInheritances.txt".format(CURRENT_DATE_VERSION), "\t")

	print "Parsing synopses..."
	phenoMimsToGenesDict = makeMIMdict()
	geneSets = defaultdict(set)

	for _, synopsis in allSynopses.iteritems():
		for keyWord in keyWords:

			## Get gene lists + disorder information, if a certain condition is met.
			if isKeyWordInAnySymptom(keyWord.lower(), synopsis):

				associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
				# brainGenes = []
				# associatedGenes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)

				geneSets[keyWord] |= set(associatedGenes)

	# Write out each gene set.
	for geneSetName, geneSet in geneSets.items():
		with open(outputFolder + geneSetName + ".tsv", "w") as out:
			out.write("\n".join(sorted(geneSet)))

	return

# Return a dict of keyterms to genes.
def getGeneSetsFromFolder(path):

	geneSets = {} # defaultdict(set)

	filenames = os.listdir(GTRfilesPath)

	for filename in filenames:

		genes = {gene.strip() for gene in readlines(filename)}
		geneSets[filename] = genes

	return geneSets


def compareSets():

	GTRgeneSets = getGTRgeneSets(GTRfilesPath)
	ourGeneSets = getGTRgeneSets(GTRfilesPath)
	assert(len(GTRgeneSets) == len(ourGeneSets))
	keywords = GTRgeneSets.keys()

	with open("Results/FirstPassComparison.tsv", "w"):
		for keyword in keyWords:

			GTRgenes = GTRgeneSets[keyword]
			ourGenes = ourGeneSets[keyword]






if __name__ == '__main__':


	# OMIM-scraping. Do this every few months.
	# scrapeOmim()
	# createOneBigFile()
	# makeInheritanceFile()
	# makeInheritanceAndSystemsFile()
	# getGeneListsForEachCategory()
	# getSymptomsFromAllCategories()

	# getMultipleSpecifiedListsOfGenes()
	# getKeyTerms()

	keyWords = getKeyTerms()
	outputFolder = PROJECT + "Results/GeneSets/Ours/"
	createOurGeneSets(keyWords, outputFolder)
	
	
	