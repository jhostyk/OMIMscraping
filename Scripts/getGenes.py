# -*- coding: utf-8 -*-
### The working functions

### Joseph Hostyk

### 2018-4-9

from helper import *


# New  version: get genes from specified categories.
# (Original version in helper.py - does it by disorder.)
def getMultipleSpecifiedListsOfGenes():

	allSynopses = getAllSynopses()
	# IDD = makeSetFromFile(PROJECT + "idd/broadIddDisordersFromOMIM.txt")
	# brainExpressedGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")
	# inheritanceDictionary = makeDictFromTwoColumns(PROJECT + "allGenes/{}/allGenesInheritances.txt".format(CURRENT_DATE_VERSION), "\t")
	iddGenes = set()

	print "Parsing synopses..."
	phenoMimsToGenesDict = makeMIMdict()
	# geneInfo = defaultdict(set)
	# keyWord = "acular degen"

	geneInfo = defaultdict(lambda: defaultdict(set))
	# Josh Cook:
	# keyWords = ["glucose intolerance", "diabetes", "insulin resistance", "acanthosis nigricans", "acromegaly", "growth hormone excess", "hypercortisolism", "cortisol excess", "cushing", "glucagonoma", "hyperinsulinemia", "hyperinsulinism", "maturity onset diabetes of the young", "neonatal diabetes", "insulitis", "beta cell", "obesity", "obese", "overweight", "pancreatitis"]
	# Kate:
	# keyWords = ["stillbirth", "stillborn", "embryonic lethality", "perinatal lethality", "neonatal lethality", "fetal lethality", "fetal demise", "embryonic demise", "perinatal demise", "neonatal demise", "perinatal death", "embryonic death", "fetal death", "neonatal death", "onset in infancy", "neonatal onset", "onset at birth", "onset in utero", "sudden death within first days of life", "sudden infant death", "sudden  cardiac death"]
	# keyWords = ["stillbirth", "stillborn", "embryonic lethality", "perinatal lethality", "neonatal lethality", "fetal lethality fetal demise", "embryonic demise", "perinatal demise", "neonatal demise", "perinatal death", "embryonic death", "fetal death", "neonatal death"]
	
	keyWords = ["inflammation", "inflammatory", "inflam"]


	# wantedCategories = ["cardiovascularHeart", "cardiovascular", "cardiovascularVascular", "respiratoryAirways", "respiratoryLarynx", "respiratoryLung", "respiratoryNasopharynx", "respiratory"]

	bronchiolitisGenes = set()
	for _, synopsis in allSynopses.iteritems():
		for keyWord in keyWords:


		## Get gene lists across a set of categories:
	# 	for category in wantedCategories:
	# 		if category in synopsis:
	# 			bronchiolitisGenes |= set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict))
	# print "There are {} bronchiolitis genes.".format(len(bronchiolitisGenes))

			## Get gene lists + disorder information, if a certain condition is met.
			if isKeyWordInAnySymptom(keyWord, synopsis):
			# if wantedCategory in synopsis:
				associatedGenes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict))
				brainGenes = []
				associatedGenes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)

				for gene in associatedGenes:

					systems = " | ".join(getSystemsAndSymptoms(synopsis, []))
					# geneInfo[gene].add((synopsis["preferredTitle"], actualInheritance, systems))

					geneInfo[gene]["phenotypes"].add(synopsis["preferredTitle"])
					geneInfo[gene]["flaggedKeywords"].add(keyWord)

	# Individual files:
	# geneSets = [geneInfo]
	# geneSetNames = [keyWord]
	# for geneSetInfo, geneSetName in zip(geneSets, geneSetNames):
	# for geneSetName, geneSetInfo in geneInfo.items():
	# 	with open(PROJECT + "FilesForJosh/" + geneSetName + "Genes.tsv", "w") as out:
	# 		out.write("2019-9-23: Scraped for {}.\n".format(geneSetName))
	# 		# out.write("Gene\tAssociated Disorder(s)\t\n")
	# 		# for gene in sorted(geneSetInfo.keys()):
	# 		# 	for associatedPhenotype in geneSetInfo[gene]:
	# 		# 		out.write("{}\t{}\n".format(gene, " | ".join(geneSetInfo[gene])))
	# 		out.write("Gene\tAssociated Disorder\tInheritance\tSystems: Symptoms\n")
	# 		for gene in sorted(geneSetInfo.keys()):
	# 			for (associatedPhenotype, inheritance, systems) in geneSetInfo[gene]:
	# 				out.write("{}\t{}\t{}\t{}\n".format(gene, associatedPhenotype, inheritance, systems))

	# Do it in one file:
	with open(PROJECT + "Inflammation.tsv", "w") as out:
		out.write("2019-9-23: Scraped for {}.\n".format(" | ".join(keyWords)))
		out.write("Gene\tAssociated Disorder(s)\tFlagged Keyword(s)\n")

		for gene, geneInfo in geneInfo.items():

			phenotypes = geneInfo["phenotypes"]
			flaggedKeywords = geneInfo["flaggedKeywords"]
			out.write("{}\t{}\t{}\n".format(gene, " | ".join(phenotypes), " | ".join(flaggedKeywords)))
	return



###################################################################################################

### Gene clustering project.

###################################################################################################


# Starting off imagining just 2 synopses in the list:
def getSymptomSimilarity(allSynopses, synopsisIDlist):

	id1 = synopsisIDlist[0]
	id2 = synopsisIDlist[1]

	shared = 0.0
	total = 0.0
	dontWant = ["inheritanceExists", "oldFormatExists"]
	for key in allSynopses[id1]:
		if key not in dontWant and isinstance(allSynopses[id1][key], bool):
			total += 1
			if allSynopses[id1][key] == allSynopses[id2][key]:
				shared += 1
	return shared/total

def getComplexGenesSymptomsSimilarities():

	allSynopses = getAllSynopses()
	with open(PROJECT + "allGenes/genesToDisorders.txt", "r") as genes:
		header = genes.readline().strip().split("\t")
		print header
		nameIndex = header.index("Gene")
		idIndex = header.index("Synopsis ID(s)")
		numDisordersIndex = header.index("Num Disorders")
		for gene in genes:
			gene = gene.strip().split("\t")
			if gene[numDisordersIndex] == "2":
				ids = gene[idIndex].split(" | ")
				percentSimilar = getSymptomSimilarity(allSynopses, ids)
				print "{} has 2 symptoms with {} similarity.".format(gene[nameIndex], percentSimilar)

# 2018-6-5
def getLiverGenes():

	allSynopses = getAllSynopses()

	print "Parsing synopses..."
	phenoMimsToGenesDict = makeMIMdict()

	geneDictTemplate = lambda:{"Associated Disease(s)": set(), "Flagged Categories": set()}
	liverGenes = defaultdict(geneDictTemplate)

	categories = set(["abdomen", "abdomenExternalFeatures", "abdomenLiver", "abdomenPancreas", "abdomenBiliaryTract", "abdomenSpleen", "abdomenGastrointestinal"])
	oldCategories = set(["GI", "Liver"])
	for _, synopsis in allSynopses.iteritems():
		genes = getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict)
		for category in categories:
			if category in synopsis:
				for gene in genes:
					liverGenes[gene]["Associated Disease(s)"].add(synopsis["preferredTitle"])
					liverGenes[gene]["Flagged Categories"].add(category)
		if "oldFormat" in synopsis:
			for category in oldCategories:
				if category in synopsis["oldFormat"]:
					for gene in genes:
						liverGenes[gene]["Associated Disease(s)"].add(synopsis["preferredTitle"])
						liverGenes[gene]["Flagged Categories"].add(category)
	print liverGenes

	with open(PROJECT + "liver/liverGenes2018-6-5.txt", "w") as out:
		out.write("Genes scraped from OMIM clinical synopses that contained symptoms for the following categories:\n")
		out.write("abdomen, abdomenExternalFeatures, abdomenLiver, abdomenPancreas, abdomenBiliaryTract, abdomenSpleen, abdomenGastrointestinal.\n")
		out.write("oldFormat categories: GI, Liver.\n")
		out.write("Gene\tAssociated Disease(s)\tFlagged Categories\n")
		for gene in liverGenes:
			diseases = " | ".join(liverGenes[gene]["Associated Disease(s)"])
			categories = " | ".join(liverGenes[gene]["Flagged Categories"])
			out.write("{}\t{}\t{}\n".format(gene, diseases, categories))





if __name__ == '__main__':


	# OMIM-scraping. Do this every few months.
	# scrapeOmim()
	# createOneBigFile()
	# makeInheritanceFile()
	# makeInheritanceAndSystemsFile()
	# getGeneListsForEachCategory()
	# getSymptomsFromAllCategories()

	getMultipleSpecifiedListsOfGenes()
	
	# getQualifyingDisordersFromFiles()
	# switchDisordersToGenes()
	# compareOldAndNewGenes()


	# neurologicExistsFile = PROJECT + "neurologicExistsAssociations.txt"
	# neurologicManifestationsExistFile = PROJECT + "neurologicManifestationsExistAssociations.txt"
	# allGenesFile = PROJECT + "allGenesAssociations.txt"
	# analyzeQualifyingResults(neurologicManifestationsExistFile)

	# countNeuroGenes(allGenesFile)

	# getSymptomsFromAllCategories()	


	# graphSetData()

	# graphPLIscores()
	# findHowManyUniqueDisordersByPLI()

	# makeInheritanceAndSystemsFile()

	# getComplexGenesSymptomsSimilarities()

	# makeNickStongOMIMfile()


	# 2019-4-11: Trying to nail down the gene-recognition from strings.
	# testGeneScraper()

	# makGundiOMIMfile()

