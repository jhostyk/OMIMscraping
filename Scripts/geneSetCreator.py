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
	keyTerms = [filename.replace(".tsv", "") for filename in filenames if ".tsv" in filename]

	print ("We have {} keyterms.".format(len(keyTerms)))
	return keyTerms

def createOurGeneSets(keyWords, outputFolder):

	allSynopses = getAllSynopses()
	officialGeneList = getOfficialGenes()
	# inheritanceDictionary = makeDictFromTwoColumns(PROJECT + "allGenes/{}/allGenesInheritances.txt".format(CURRENT_DATE_VERSION), "\t")

	print ("Parsing synopses...")
	phenoMimsToGenesDict = makeMIMdict()
	geneSets = defaultdict(set)

	for _, synopsis in allSynopses.items():
		for keyWord in keyWords:
			# if keyWord == "Alzheimer's Disease":
			# 	print ("hm")
			# 	raise

			## Get gene lists + disorder information, if a certain condition is met.
			if iskeyWordInAnySymptom(keyWord.lower(), synopsis):

				associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
				# brainGenes = []
				# associatedGenes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)

				geneSets[keyWord] |= set(associatedGenes)

	# Write out each gene set.
	# print ("Heart attack:", geneSets["Heart attack"])
	for keyWord in keyWords:
		genes = geneSets[keyWord]
		with open(outputFolder + keyWord + ".tsv", "w") as out:
			out.write("\n".join(sorted(genes)))

	return

# Return a dict of keyterms to genes.
def getGeneSetsFromFolder(path):

	geneSets = {} # defaultdict(set)

	filenames = os.listdir(path)

	for filename in filenames:

		if ".tsv" in filename:

			with open(path + filename, "r") as file:
				genes = {gene.strip() for gene in file}
				geneSets[filename] = genes

	return geneSets


def compareSets():

	GTRfilesPath = PROJECT + "Results/GeneSets/GTRgeneLists/"
	ourFilesPath = PROJECT + "Results/GeneSets/Ours/"

	GTRgeneSets = getGeneSetsFromFolder(GTRfilesPath)
	ourGeneSets = getGeneSetsFromFolder(ourFilesPath)

	# print(set(GTRgeneSets.keys()).difference(set(ourGeneSets.keys())))
	# print(set(ourGeneSets.keys()).difference(set(GTRgeneSets.keys())))

	assert(len(GTRgeneSets) == len(ourGeneSets))
	# raise
	keyWords = GTRgeneSets.keys()

	with open("Results/FirstPassComparison.tsv", "w") as out:

		out.write("Set Name\tNum shared\tNumber just in GTR\tJust GTR\tNumber just in ours\tJust ours\n")
		for keyWord in sorted(keyWords):

			GTRgenes = GTRgeneSets[keyWord]
			ourGenes = ourGeneSets[keyWord]

			shared = ourGenes.intersection(GTRgenes)
			justGTR = GTRgenes.difference(ourGenes)
			justOurs = ourGenes.difference(GTRgenes)

			out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(keyWord, len(shared), len(justGTR), " | ".join(justGTR), len(justOurs), " | ".join(justOurs)))


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

	# keyWords = getKeyTerms()

	# outputFolder = PROJECT + "Results/GeneSets/Ours/"
	# createOurGeneSets(keyWords, outputFolder)

	compareSets()

	
	
	