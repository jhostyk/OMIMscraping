# -*- coding: utf-8 -*-
### The working functions

### Joseph Hostyk

### 2018-4-9

from helper import *


###################################################################################################

### First steps: compare our results to GTR's.

###################################################################################################

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

			## Get gene lists + disorder information, if a certain condition is met.
			if isKeyWordInAnySymptom(keyWord.lower(), synopsis):

				associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
				# brainGenes = []
				# associatedGenes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)

				geneSets[keyWord] |= set(associatedGenes)

	# Write out each gene set.
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

	assert(len(GTRgeneSets) == len(ourGeneSets))
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



###################################################################################################

### Next: filtering categories.

###################################################################################################


def getGeneSet(keyWord):


	 question = [inquirer.List("Keyword", message = "Please enter your keyword",
					choices = ["No, these are great!"] + list(categoriesToGenes.keys()),
				),
	]
	answer = inquirer.prompt(question)
	keyWord = answer["Keyword"]



	### TODO:
	synonyms = getSynonyms(keyword)

	allSynopses = getAllSynopses()
	officialGeneList = getOfficialGenes()

	print ("Parsing synopses...")
	phenoMimsToGenesDict = makeMIMdict()
	foundGenes = set()
	categoriesToGenes = defaultdict(set)

	dontWant = ["molecularBasis", "miscellaneous", "inheritance"]



	### Get the info
	for _, synopsis in allSynopses.items():

		## Get gene lists + disorder information, if a certain condition is met.
		if isKeyWordInAnySymptom(keyWord.lower(), synopsis):

			currentCategories = set()

			associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
			foundGenes |= associatedGenes
			for aspect in synopsis:
				if isinstance(synopsis[aspect], str) and "{" in synopsis[aspect] and aspect not in dontWant:
					currentCategories.add(aspect)
			for gene in associatedGenes:
				for category in currentCategories:

					categoriesToGenes[category].add(gene)

	### Let the user remove categories:
	categoryToExclude = None
	excludedCategories = set()
	excludedGenes = set()
	while categoryToExclude !=  "No, these are great!":

		genesToExclude = categoriesToGenes[categoryToExclude]

		for geneToExclude in genesToExclude:

			excludedGenes.add(geneToExclude)
			excludedCategories.add(categoryToExclude)
			foundGenes.remove(geneToExclude)
			if categoryToExclude:
				excludedCategories.add(categoryToExclude)
			for category, genes in categoriesToGenes.items():

				categoriesToGenes[category] = {gene for gene in genes if gene != geneToExclude}
		### If we removed the last gene from a certain category, we get rid of it:
		categoriesToGenes = {category: genes for category, genes in categoriesToGenes.items() if len(genes) > 0}

		question = [
	 	inquirer.List("Exclude", message = "There are currently {} genes. Would you like to exclude any categories?".format(len(foundGenes)),
	 				choices = ["No, these are great!"] + list(categoriesToGenes.keys()),
					),
		]
		answer = inquirer.prompt(question)
		categoryToExclude = answer["Exclude"]

	print ("{} categories were excluded, removing {} total genes.".format(len(excludedCategories), len(excludedGenes)))
	print ("The remaining categories are:")
	print (" | ".join(sorted(categoriesToGenes.keys())))
	print ("Final gene list:")
	print (sorted(foundGenes))

				

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

	getGeneSet("acanthosis")

	
	
	