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


def getGeneSet():

	conceptsToTheirCode, codesToSynonyms = getAllSynonyms()

	allSynopses = getAllSynopses()
	officialGeneList = getOfficialGenes()

	print ("Parsing synopses...")
	phenoMimsToGenesDict = makeMIMdict()

	responses = set()
	response = None

	while response != "done":

		question = [inquirer.Text("Keyword", message = "Please enter your search term, or 'done'"
					),
		]
		answer = inquirer.prompt(question)
		response = answer["Keyword"]

		if response != "done":
			responses.add(response)
	print("using the keyWords:", responses)

	keyWords = set(responses)
	for response in responses:
		### If a term doesn't have a synonym, it won't be in conceptsToTheirCode
		if response.lower() in conceptsToTheirCode:
			code = conceptsToTheirCode[response.lower()]
			keyWords |= set(codesToSynonyms[code])
			print ("For {}, found the synonyms {}.".format(response, codesToSynonyms[code]))
	
	# foundGenes = set()
	foundGenesToKeyWords = defaultdict(set)
	categoriesToGenes = defaultdict(set)

	dontWant = ["molecularBasis", "miscellaneous", "inheritance"]



	### Get the info
	for _, synopsis in allSynopses.items():
		for keyWord in keyWords:

			## Get gene lists + disorder information, if a certain condition is met.
			if isKeyWordInAnySymptom(keyWord.lower(), synopsis):

				currentCategories = set()

				associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
				# foundGenes |= associatedGenes
				for gene in associatedGenes:
					foundGenesToKeyWords[gene].add(keyWord)
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
			# foundGenes.remove(geneToExclude)
			del foundGenesToKeyWords[geneToExclude]
			if categoryToExclude:
				excludedCategories.add(categoryToExclude)
			for category, genes in categoriesToGenes.items():

				categoriesToGenes[category] = {gene for gene in genes if gene != geneToExclude}
		### If we removed the last gene from a certain category, we get rid of it:
		categoriesToGenes = {category: genes for category, genes in categoriesToGenes.items() if len(genes) > 0}

		question = [
	 	inquirer.List("Exclude", message = "There are currently {} genes. Would you like to exclude any categories?".format(len(foundGenesToKeyWords)),
	 				choices = ["No, these are great!"] + list(categoriesToGenes.keys()),
					),
		]
		answer = inquirer.prompt(question)
		categoryToExclude = answer["Exclude"]

	# print ("{} categories were excluded, removing {} total genes.".format(len(excludedCategories), len(excludedGenes)))
	# print ("The remaining categories are:")
	# print (" | ".join(sorted(categoriesToGenes.keys())))
	print ("Final gene list:")
	print (sorted(foundGenesToKeyWords))

	with open(PROJECT + "test.tsv", "w") as out:

		out.write("Gene\t{}\n".format("\t".join(sorted(keyWords))))
		for gene, keyWordsThatFlaggedIt in sorted(foundGenesToKeyWords.items()):

			out.write("{}\t{}\n".format(gene, "\t".join(["X" if keyWord in keyWordsThatFlaggedIt else "" for keyWord in sorted(keyWords)])))

				

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

	getGeneSet()
	# makeUniqueCodes()

	
	
	