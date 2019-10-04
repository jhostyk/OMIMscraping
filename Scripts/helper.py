# -*- coding: utf-8 -*-
### Helper functions and older versions of functions.

### Joseph Hostyk

### 2018-1-12

PROJECT = "/nfs/labs/goldstein/jh3958/NeuroGenesPaper/"
CURRENT_DATE_VERSION = "2019-08-02"
# import datetime
# CURRENT_DATE_VERSION = datetime.datetime.now().strftime("%Y-%m-%d")


from OMIMapiKey import *


import os
import sys
import time
import urllib2
import imp
import re
from difflib import SequenceMatcher
from collections import defaultdict

###################################################################################################

### OMIM Scraping

###################################################################################################

# Make a call and save 20 OMIM results.
def saveOmimResults(startIndex):

	endIndex = str(startIndex + 20)
	startIndex = str(startIndex)
	pythonUrl = "https://api.omim.org/api/clinicalSynopsis/search?search=*&filter=&fields=&start={}&limit={}&sort=&operator=&include=clinicalSynopsis&include=existFlags&format=python&apiKey={}".format(startIndex, endIndex, API_KEY)
	response = urllib2.urlopen(pythonUrl)
	webContent = response.read()

	title = PROJECT + "omimResults/omimResults{}to{}.py".format(startIndex, endIndex)
	f = open(title, 'w')
	f.write("omim = ")
	f.write(webContent)
	f.close

# Go through all OMIM results.
def scrapeOmim():

	totalResults = 7368
	resultsLimit = 20
	numSearches = totalResults/resultsLimit
	for i in range(numSearches+1):
		sys.stdout.write("Scraping webpage {} of {}\r".format(i, numSearches))
		sys.stdout.flush()
		# OMIM only gives 20 results at a time.
		saveOmimResults(i * 20)
		time.sleep(3)

# Had been working with a ton of files. Instead, pull out the important info,
# concat, and write to one big file.
def createOneBigFile():

	sys.dont_write_bytecode = True
	path = PROJECT + "omimResults/"
	outputFileName = path + "{}_all.py".format(CURRENT_DATE_VERSION)
	files = [file for file in os.listdir(path) if "omimResults" in file]
	numSynopses = 0
	with open(outputFileName, "w") as outfile:
		outfile.write("omim = {\n")
		for file in files:
			with open(path + file, "r") as infile:
				# First 20 lines are metadata
				for _ in range(19):
					infile.readline()
				
				for line in infile:
					if line == "{'clinicalSynopsis': { \n":
						outfile.write("\"clinicalSynopsis{}\": {{\n".format(numSynopses))
						numSynopses += 1
						continue

					# This bit marks the end of each synopsis.
					skipStrings = ["} },\n", "} \n", "} ] \n", "} \n", "} }"]
					if any(line == skip for skip in skipStrings):
						continue
					if "'matches':" in line:
						outfile.write(line.strip())
						outfile.write("}, \n")
						continue
					outfile.write(line)
		outfile.write("}")		

###################################################################################################

### Helper functions

###################################################################################################

def getAllSynopses():

	print "Loading data..."
	sys.dont_write_bytecode = True
	allSynopses = imp.load_source('', PROJECT + "omimResults/{}_all.py".format(CURRENT_DATE_VERSION)).omim

	return allSynopses


def makeSetFromFile(filename):

	stuff = set()
	with open(filename, "r") as file:
		file.readline() # skip the header
		for line in file:
			stuff.add(line.strip())
	return stuff

def makeDictFromTwoColumns(filename, delimiter):

	stuff = dict()
	with open(filename, "r") as file:
		file.readline() # skip the header
		for line in file:
			line = line.strip().split(delimiter)
			stuff[line[0]] = line[1]
	return stuff


# Initial step just to figure out the structure of the OMIM-scraped Python files.
def dictParser(dic, numTabs):

	for key, value in dic.iteritems():
		print "{}key: {}".format("\t"*numTabs, key)
		if isinstance(value, dict):
			dictParser(value, numTabs + 1)
		else:
			print "{}Value: {}".format("\t"*numTabs, value)

# dictParser(omim["omim"]["searchResponse"]["clinicalSynopsisList"][0], numTabs)#["neurologicExists"]

# Lets us find associated genes more easily, later.
def makeMIMdict():

	phenoMimsToGenesDict = {}
	
	with open(PROJECT + "ReferenceFiles/genemap2.txt", "r") as mims:
		geneMimIndex = 5
		officialGeneNameIndex = 8
		allGeneNamesIndex = 6
		phenoIndex = 12
		for line in mims:
			# Skip the commented lines
			if line[0] == "#":
				continue
			line = line.split("\t")
			geneMim = line[geneMimIndex]
			gene = line[officialGeneNameIndex]
			if len(gene) == 0:
				# In the "gene symbols" column, first name is usually the right one if nothing 
				# is listed in the official column.
				gene = line[allGeneNamesIndex].split(", ")[0]
			pheno = line[phenoIndex]



			# Annoying, but the only way to pull out the phenotype MIM (not the gene MIM).
			# We make a set because there are sometimes duplicate MIMs.
			phenoMims = set(re.findall(r"(\d+) \(\d\)", pheno))
			if len(gene) == 0 and len(phenoMims) == 0:
				continue
				# print "phenomim {} has a jank gene.".format(phenoMims)
				# if len(phenoMims) == 0:
				# 	print "Gene: {}, geneMim: {}, pheno: {}, regex phenoMim: {}".format(gene, geneMim, pheno, phenoMims)

			for phenoMim in phenoMims:
				if phenoMim not in phenoMimsToGenesDict:
					phenoMimsToGenesDict[phenoMim] = []
				phenoMimsToGenesDict[phenoMim] += [(gene, geneMim)]

	return phenoMimsToGenesDict

# Original version. Doesn't miss anything, but adds in a lot of junk.
# Pull out the gene name.
# def getGenesFromMolecularBasisString(molecularBasis):
# 	# Original:
# 	pattern = re.compile(r"\((\w*)\, \{")
# 	genes = []
# 	split = molecularBasis.split("gene")
# 	for section in split:
# 		result = pattern.findall(section)
# 		if result:
# 			genes += result
# 	return genes
# Isn't getting genes from these strings:
# Caused by mutation in the complement component 3 gene (C3
# the SLC3A1 ({104614}), PREPL ({609557}), PPM1B ({603770}), and C2orf34 ({609559}) genes
# Caused by mutation in the MMACHC gene (MMACHC gene {609831.0001})
# Caused by mutation in the PKHD1 gene ({606702.0001})
# Caused by mutation in the glypican 3 gene (GPC3, 300037.0001) {UMLS C1839289}

# New version: has two methods.
def getGenesFromMolecularBasisString(molecularBasis):
	
	potentialGenes = set()

	### Check that all caps and numbers:
	molecularBasisWords = molecularBasis.replace("(", "").replace(",", "").split(" ")
	GOOD_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
	for word in molecularBasisWords:
		def isGood(string):
			for letter in string:
				if letter not in GOOD_CHARS:
					return False
			return True
		def isJustNumbers(string):
			try:
				int(string)
				return True # If it can be turned into an int, there are no letters.
			except ValueError as e:
				return False
		if isGood(word) and not isJustNumbers(word):
			potentialGenes.add(word)
	# Misses things like "C2orf34".

	### Check based on string position.

	molecularBasisWords = molecularBasis.split(" ")
	for word in molecularBasisWords:
		if word[0] == "(" and word[-1] == ",":
			potentialGenes.add(word.replace("(","").replace(",", ""))
	
	return list(potentialGenes)

	# Won't catch "[] gene", e.g. SHANK3.
	# Search for "({" when looking for missed genes.
	
def testGeneScraper():

	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()

	print "Parsing synopses..."
	genesInheritance = {}
	allMaybeGenes = set()
	for _, synopsis in allSynopses.iteritems():
		if "molecularBasis" in synopsis:
			molecularBasis = synopsis["molecularBasis"]
		else:
			molecularBasis = "No molecular basis listed in OMIM's clinical synopsis."
		allMaybeGenes |= set(TESTgetGenesFromMolecularBasisString(molecularBasis))
	with open(PROJECT + "testGenes.txt", "w") as out:
		out.write("\n".join(allMaybeGenes))

def makeGeneListPretty(array):

	stringArray = str(array)
	clean = stringArray.replace("[\'", "").replace("\', \'", ", ").replace("\']", "").replace("\"", "")

	return clean

# def getGeneListFromMimAndMolecularBasis(phenoMim, molecularBasis):
def getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict):

	
	if "molecularBasis" in synopsis:
		molecularBasis = synopsis["molecularBasis"]
	else:
		molecularBasis = "No molecular basis listed in OMIM's clinical synopsis."

	phenoMim = str(synopsis["mimNumber"])

	try:
		mapGenes = [gene for (gene, _) in phenoMimsToGenesDict[phenoMim]]
	except KeyError as e:
		mapGenes = "No associated gene found in OMIM's genemap2.txt."

	molecularBasisScrapedGenes = getGenesFromMolecularBasisString(molecularBasis)


	prettyMapGenes = makeGeneListPretty(mapGenes)
	prettyMolecularBasisGenes = makeGeneListPretty(molecularBasisScrapedGenes)


	genes = set()
	if len(prettyMapGenes) > 2 and "No associated gene" not in prettyMapGenes:
		genes |= (set(prettyMapGenes.split(", ")))
	if len(prettyMolecularBasisGenes) > 2:
		genes |= (set(prettyMolecularBasisGenes.split(", ")))
	if len(prettyMolecularBasisGenes) < 3 and  "No molecular basis" not in molecularBasis:# and "Caused by mutation" not in molecularBasis:
		genes.add(molecularBasis)

	if len(genes) == 0:
		genes.add("No known molecular basis.")
	genes = makeGeneListPretty(list(genes)).split(", ")

	return genes


def isKeyWordInAnySymptom(keyword, synopsis):

	for key in synopsis:
		if isinstance(synopsis[key], basestring) and keyword in synopsis[key].lower():
			return True
	if "oldFormat" in synopsis:
		return isKeyWordInAnySymptom(keyword, synopsis["oldFormat"])
	return False

def getSystemsAndSymptoms(synopsis, systemsList):

	dontWant = ["prefix", "matches", "preferredTitle", "inheritance", "Inheritance"]
	for key in synopsis:
		if isinstance(synopsis[key], basestring) and key not in dontWant:
			systemsList.append("{}: {}".format(key, synopsis[key].replace("\n", " "	)))
	if "oldFormat" in synopsis:
		return getSystemsAndSymptoms(synopsis["oldFormat"], systemsList)
	return systemsList

	
###################################################################################################

### The actual important functions.

###################################################################################################

def getQualifyingDisordersFromDict(omim):
	qualifyingPhenoMims = []
	# categories = ["neurologicExists"]
	categories = ["neurologicBehavioralPsychiatricManifestationsExists"]
	# categories = ["genitourinaryKidneysExists"]
	disorders = omim["omim"]["searchResponse"]["clinicalSynopsisList"]
	
	for disorder in disorders:
		title = disorder["clinicalSynopsis"]["preferredTitle"]
		mim = str(disorder["clinicalSynopsis"]["mimNumber"])
		# Check if the disorder matches our criteria
		for category in categories:
			if disorder["clinicalSynopsis"][category]:
				if "molecularBasis" in disorder["clinicalSynopsis"]:
					molecularBasis = disorder["clinicalSynopsis"]["molecularBasis"]
				else:
					molecularBasis = "No molecular basis listed in OMIM's clinical synopsis."
				qualifyingPhenoMims += [(mim, title, molecularBasis)]
	return qualifyingPhenoMims


def getQualifyingDisordersFromFiles():

	sys.dont_write_bytecode = True
	path = PROJECT + "omimResults/"
	files = os.listdir(path)
	qualifyingPhenoMims = []
	for file in files:
		foo = imp.load_source('', path + file)
		qualifyingPhenoMims += getQualifyingDisordersFromDict(foo.omim)

	# with open(PROJECT + "neurologicExistsAll.txt", "w") as out:
	# with open(PROJECT + "genitourinaryKidneysAll.txt", "w") as out:
	with open(PROJECT + "Disorders and Molecular Bases.txt", "w") as out:
		out.write("Disorder\tMolecular Basis\n")
		for phenoMim, title, molecularBasis in qualifyingPhenoMims:
			genes = getGeneListFromMimAndMolecularBasis(phenoMim, molecularBasis)
			out.write("{}\t{}\n".format(title, genes))#, prettyMolecularBasisGenes, molecularBasis.replace("\n", " ")))
	return


def switchDisordersToGenes():

	with open(PROJECT + "Disorders and Associated Genes.txt", "r") as inputFile:
		with open(PROJECT + "Genes and Associated Disorders.txt", "w") as out:
			inputFile.readline()
			out.write("Gene\tAssociated Disorder\n")
			disorderIndex = 0
			genesIndex = 1

			for line in inputFile:
				line = line.split("\t")
				if "No known molecular basis." not in line[genesIndex]:
					genes =  makeGeneListPretty(line[genesIndex].strip()).split(", ")
					# print genes
					for gene in genes:
						out.write("{}\t{}\n".format(gene, line[disorderIndex]))
				
#################
### Written after put all the synopses in one file, so not as clunky as code above.
#################

###################################################################################################

### Data exploration: Neuro.

###################################################################################################


def getManifestationsFromString(string):

	splitString = string.split("\n")
	splitString = ["FLAG" + s for s in splitString]
	pattern = re.compile(r"FLAG([\w ]*) \{")
	for s in splitString:
		manifestations = pattern.findall(s)
	return manifestations


def getManifestationStringFromSynopsis(synopsis, category, oldFormatCategory):

	manifestations = set()
	if category in synopsis:
		manifestationsString = synopsis[category]
		manifestations |= set(getManifestationsFromString(manifestationsString))
	if "oldFormat" in synopsis:
		manifestationsString = synopsis["oldFormat"][oldFormatCategory]
		manifestations |= set(getManifestationsFromString(manifestationsString))
	
	return manifestations


def getManifestations():

	allSynopses = getAllSynopses()	
	
	# manifestations = set()
	categories = [('genitourinaryKidneys', "Renal"), ('headAndNeckEyes', "Eyes"), ('abdomenLiver', "Liver")]
	kidney = set()
	eyes = set()
	liver = set()

	for synopsisID, synopsis in allSynopses.iteritems():
		kidney |= getManifestationStringFromSynopsis(synopsis, category = "genitourinaryKidneys", oldFormatCategory = "Renal")
		eyes |= getManifestationStringFromSynopsis(synopsis, category = "headAndNeckEyes", oldFormatCategory = "Eyes")
		liver |= getManifestationStringFromSynopsis(synopsis, category = "abdomenLiver", oldFormatCategory = "Liver")

	for (category, oldFormat), symptomSet in zip(categories, [kidney, eyes, liver]):
		with open(PROJECT + "GeneSetsForKate/" + category + "Symptoms.txt", "w") as out:
			out.write("Symptoms found for {} (old format category: '{}'), in OMIM's clinical synopses.\n".format(category, oldFormat))
			out.write("\n".join(symptomSet))



# 2018-6-4: To match up symptoms better, get all the oldFormat symptoms.
def getAllOldFormatSymptoms():

	allSynopses = getAllSynopses()
	oldSynopses = set()


def findOldFormatMappings():

	oldFormatToNew = defaultdict(set)
	allSynopses = getAllSynopses()
	for _, synopsis in allSynopses.iteritems():
		if "oldFormat" in synopsis:
			allComorbidities = [category for category in synopsis if synopsis[category] == True]
			for oldCategory in synopsis["oldFormat"]:
				oldFormatToNew[oldCategory] |= set(allComorbidities)
	print oldFormatToNew

# 2018-6-4: Write h
def getSymptomsFromAllCategories():

	dontWant = ["prefix", "matches", "preferredTitle", "inheritance", "Inheritance"]

	categoriesToSymptoms = defaultdict(set)
	allSynopses = getAllSynopses()
	for _, synopsis in allSynopses.iteritems():
		# if "oldFormat" in synopsis:
			# oldFormatDict = synopsis["oldFormat"]
			# for category in oldFormatDict:
			# 	categoriesToSymptoms[category] |= getManifestationStringFromSynopsis(oldFormatDict, category, oldFormatCategory=None)
		for category in synopsis:
			if isinstance(synopsis[category], str) and category not in dontWant:
				categoriesToSymptoms[category] |= getManifestationStringFromSynopsis(synopsis, category, oldFormatCategory=None)
		# print categoriesToSymptoms
		# raise
	# print categoriesToSymptoms
	print "Writing results:"
	for category in categoriesToSymptoms:
		with open(PROJECT + "ClinicalSynopsisCategories/{date}/{cat}.txt".format(date = CURRENT_DATE_VERSION, cat = category), "w") as out:
			out.write("{date}. All symptoms listed in category {cat}:\n".format(date = CURRENT_DATE_VERSION, cat = category))
			out.write("\n".join(categoriesToSymptoms[category]))


###################################################################################################

### Creating files.

###################################################################################################

def makeIDDdisordersFile():

	iddSymptoms = {"neurologicCentralNervousSystem": ["Intellectual disability", "Mental retardation", "Delayed development", "Global developmental delay", "Globally delayed development", "Autism", "Mild cognitive impairment", "Cognitive deficits", "Severe mental retardation", "Autistic behavior", "Autistic features", "Developmental delay", "Neurodevelopmental regression", "Mild mental retardation", "Delayed cognitive development", "Mental impairment", "Variable mental retardation ranging from severe neurodegeneration to mild mental retardation", "Mental impairment may develop with repeated acute episodes", "Mental retardation has been reported", "Developmental disability", "Normal to mild mental retardation", "Profound mental retardation and hypotonia in survivors", "Mental retardation due to repeated episodes of hypoglycemia", "Mental deterioration", "Mental retardation can occur in patients with repeated episodes of dehydration", "Global developmental delay if untreated", "Occasional mental retardation", "Autistic behaviors", "Cognitive delay", "Cognitive impairment"], "neurologicBehavioralPsychiatricManifestations": ["Autism", "Autistic"] }

	allSynopses = getAllSynopses()

	print "Parsing synopses..."
	idds = []
	for _, synopsis in allSynopses.iteritems():
		for symptomType in iddSymptoms:
			if symptomType in synopsis:
				for symptom in iddSymptoms[symptomType]:
					if symptom in synopsis[symptomType]:
						idds.append(synopsis["preferredTitle"])

	with open(PROJECT + "idd/broadIddDisordersFromOMIM.txt", "w") as out:
		out.write("Scraped all OMIM synopses, in 'Neurologic/central nervous system', for the following symptoms: {}.\n".format(iddSymptoms))
		for idd in idds:
			out.write(idd + "\n")

def makeBrainExpressedGenesFile():
	path = "/nfs/seqscratch11/jh3958/GTEx_Analysis_v7_eQTL/"
	files = os.listdir(path)
	qualifyingPhenoMims = []
	genesExpressedInNonBrain = set()
	brainGenes = set()
	for file in files:
		if ".egenes.txt" in file:
			print "Reading ", file
			with open(path + file, "r") as gtex:
				header = gtex.readline().split("\t")
				geneIndex = header.index("gene_name")
				for line in gtex:
					line = line.split("\t")
					if "Brain" in file:
						brainGenes.add(line[geneIndex])
					else:
						genesExpressedInNonBrain.add(line[geneIndex])
	print "Expressed in non-brain organ (might include brain ones):", len(genesExpressedInNonBrain)
	print "Expressed in brain:", len(brainGenes)
	exclusiveBrain = brainGenes.difference(genesExpressedInNonBrain)
	print "Exclusively brain:", len(exclusiveBrain)

	with open(PROJECT + "ReferenceFiles/allBrainExpressedGenesFromGTEx.txt", "w") as out:
		for gene in brainGenes:
			out.write(gene + "\n")
	with open(PROJECT + "ReferenceFiles/exclusivelyBrainExpressedGenesFromGTEx.txt", "w") as out:
		for gene in exclusiveBrain:
			out.write(gene + "\n")

# Don't know how to find genes' inheritance. Can only get per disorder.
# Want to see if for each gene, all the disorders they cause have the same type of inheritance.
def makeInheritanceFile():

	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()
	brainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")

	print "Parsing synopses..."
	genesInheritance = {}
	n = 0
	for _, synopsis in allSynopses.iteritems():
		sys.stdout.write("\rOn synopsis {}".format(n))
		sys.stdout.flush()
		n+=1
		genes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)
		
		for gene in genes:#.split(", "):
			if gene not in genesInheritance:
				genesInheritance[gene] = set()
			genesInheritance[gene].add(actualInheritance)
	with open(PROJECT + "allGenes/{}/allGenesInheritances.txt".format(CURRENT_DATE_VERSION), "w") as out:
		out.write("Gene\tInheritance\n")

		for gene in genesInheritance:
			out.write("{}\t{}\n".format(gene, ", ".join(genesInheritance[gene])))

# Same as above, but also info on the gene systems
def makeInheritanceAndSystemsFile():

	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()
	brainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")

	print "Parsing synopses..."
	genesData = {}
	dontWant = ["molecularBasis", "miscellaneous", "inheritance"]
	for _, synopsis in allSynopses.iteritems():

		genes, inheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)
		systems = set()
		for aspect in synopsis:
			if isinstance(synopsis[aspect], basestring) and "{" in synopsis[aspect] and aspect not in dontWant:
				systems.add(aspect)


		for gene in genes:
			if gene not in genesData:
				genesData[gene] = {"Phenotypes": set(), "Inheritance": set(), "Associated Systems": set()}
			genesData[gene]["Phenotypes"].add(synopsis["preferredTitle"])
			genesData[gene]["Inheritance"].add(inheritance)
			genesData[gene]["Associated Systems"] |= systems
	with open(PROJECT + "allGenes/{}/allGenesInheritancesAndSystems.txt".format(CURRENT_DATE_VERSION), "w") as out:
		out.write("Gene\tAssociated Phenotype(s)\tInheritance\tAssociated System(s)\n")
		for gene in genesData:
			out.write("{}\t{}\t{}\t{}\n".format(gene, " | ".join(genesData[gene]["Phenotypes"]), ", ".join(genesData[gene]["Inheritance"]), ", ".join(genesData[gene]["Associated Systems"])))

# 2018-10-24: Make file of Gene: Condition1 [Inheritance] | Condition2 [Inheritance]
def makeNickStongOMIMfile():

	abbreviations = {"Autosomal Dominant": "AD", "Autosomal Recessive": "AR", "X-Linked Recessive": "XLR",  "X-Linked Dominant": "XLD"}

	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()
	brainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")

	print "Parsing synopses..."
	genesData = defaultdict(set)
	for _, synopsis in allSynopses.iteritems():

		genes, inheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)

		if inheritance in abbreviations:
			inheritance = abbreviations[inheritance]

		for gene in genes:
			genesData[gene].add(synopsis["preferredTitle"] + " [{}]".format(inheritance))
	with open(PROJECT + "allGenes/{}/allGenesToAssociatedConditions.txt".format(CURRENT_DATE_VERSION), "w") as out:
		out.write("Gene\tAssociated Phenotype(s) and Inheritance\n")
		for gene in genesData:
			out.write("{}\t{}\n".format(gene, " | ".join(genesData[gene])))

# 2018-10-29: Make file where each synopsis gets it own line, per gene.
def makGundiOMIMfile():

	CURRENT_DATE_VERSION = "2019-4-12"
	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()
	brainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")

	print "Parsing synopses..."
	genesData = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
	dontWant = ["molecularBasis", "miscellaneous", "inheritance"]
	for _, synopsis in allSynopses.iteritems():

		genes, inheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict)
		systems = set()
		for aspect in synopsis:
			if isinstance(synopsis[aspect], basestring) and "{" in synopsis[aspect] and aspect not in dontWant:
				systems.add(aspect)


		for gene in genes:
			genesData[gene][synopsis["preferredTitle"]]["Inheritance"] = inheritance
			genesData[gene][synopsis["preferredTitle"]]["Systems"] = systems
	with open(PROJECT + "{}_gundiAllGenesInheritancesAndSystems.txt".format(CURRENT_DATE_VERSION), "w") as out:
		out.write("Gene\tAssociated Phenotype\tInheritance\tAssociated System(s)\n")
		for gene in genesData:
			for synopsis in genesData[gene]:
				out.write("{}\t{}\t{}\t{}\n".format(gene, synopsis, genesData[gene][synopsis]["Inheritance"], "|".join(list(genesData[gene][synopsis]["Systems"]))))


def makePLIfile():

	with open("/nfs/goldstein/software/atav_home/data/rvis/gene_score_140318.csv", "r") as genes:
		header = genes.readline().split(",")
		geneIndex = header.index("Gene")
		pliIndex = header.index("LoF-pLI[ExAC]")

		with open(PROJECT + "ReferenceFiles/pliScores.txt", "w") as out:
			out.write("Gene\tPLIscore\t(gotten from /nfs/goldstein/software/atav_home/data/rvis/gene_score_140318.csv)\n")
			for gene in genes:
				gene = gene.split(",")
				out.write("{}\t{}\n".format(gene[geneIndex], gene[pliIndex]))


def getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict):

	genes = getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict)

	# allBrainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/allBrainExpressedGenesFromGTEx.txt")
	# exclusivelyBrainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/exclusivelyBrainExpressedGenesFromGTEx.txt")
	brainExpressed = "Not expressed in brain."
	for gene in genes:
		if gene in brainGenes:
			brainExpressed = "Expressed in brain."

	if not synopsis["inheritanceExists"]:
		return (genes, "NA", brainExpressed)


	try:
		inheritance = synopsis["inheritance"]

	# Old-formatted synopses are annoying.
	except KeyError as e:
		inheritance = synopsis["oldFormat"]["Inheritance"]
		inheritance = inheritance[:inheritance.index(";")]
	if "{" in inheritance:
		inheritance = inheritance[:inheritance.index("{")]
	inheritance = inheritance.replace("?","").strip(" ")

	# actualInheritances = ["Autosomal Dominant", "Autosomal Recessive", "X-Linked Recessive", "X-Linked Dominant"]

	# actualInheritance = "Other"
	# actualInheritance = inheritance
	# for i in actualInheritances:
	# 	if SequenceMatcher(None, inheritance, i).ratio() > 0.7:
	# 		actualInheritance = i

	abbreviations = {"Autosomal Dominant": "AD", "Autosomal Recessive": "AR", "X-Linked Recessive": "XLR",  "X-Linked Dominant": "XLD"}

	actualInheritance = inheritance
	for inheritance in abbreviations:
		# The strings are not standardly formatted. Some have spaces at the end. This matching works well.
		if SequenceMatcher(None, actualInheritance, inheritance).ratio() > 0.7:
			actualInheritance = abbreviations[inheritance]

			
	return (genes, actualInheritance, brainExpressed)


# Original version: get disorders of specified categories.
def getMultipleSpecifiedListsOfDisorders():

	print "Loading data..."
	IDD = set()
	with open(PROJECT + "idd/broadIddDisordersFromOMIM.txt", "r") as idds:
		idds.readline() # The header
		for idd in idds:
			IDD.add(idd.strip())

	allSynopses = getAllSynopses()

	iddDisorders = {}
	behavioralDisorders = {}
	kidneyDisorders = {}

	print "Parsing synopses..."
	for _, synopsis in allSynopses.iteritems():

		# Take care of IDD
		notIDD = True
		for idd in IDD:
			if idd == synopsis["preferredTitle"]:
				iddDisorders[idd] = getSynopsisInfo(synopsis)
				notIDD = False
		if notIDD:
			if synopsis["neurologicBehavioralPsychiatricManifestationsExists"]:
				behavioralDisorders[synopsis["preferredTitle"]] = getSynopsisInfo(synopsis)

		if synopsis["genitourinaryKidneysExists"]:
			kidneyDisorders[synopsis["preferredTitle"]] = getSynopsisInfo(synopsis)

	disorders = [iddDisorders, behavioralDisorders, kidneyDisorders]
	disorderFileNames = ["iddData.txt", "behavioralData.txt", "kidneyData.txt"]

	print "Writing results..."
	for disorder, disorderFileName in zip(disorders, disorderFileNames):
		with open(PROJECT + disorderFileName, "w") as out:
			out.write("Disorder\tMolecular Basis\tInheritance\tBrain Expressed?\n")
			for disorder, (genes, inheritance, brainExpressed) in disorder.iteritems():
				out.write("{}\t{}\t{}\t{}\n".format(disorder, genes, inheritance, brainExpressed))
	return
	
### Making reference files, using all the scraped data.
# Get and write each gene's associated disorders.
def getAllGenesDisorders():

	allSynopses = getAllSynopses()

	phenoMimsToGenesDict = makeMIMdict()

	genesAndTheirDisorders = {}

	print "Parsing synopses..."
	for synopsisID, synopsis in allSynopses.iteritems():

		genes = getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict)
		for gene in genes:
			if gene not in genesAndTheirDisorders:
				genesAndTheirDisorders[gene] = {"titles": set(), "indices": set()}
			genesAndTheirDisorders[gene]["indices"].add(synopsisID)
			genesAndTheirDisorders[gene]["titles"].add(synopsis["preferredTitle"])

	print "Writing results..."
	
	with open(PROJECT + "allGenes/genesToDisorders.txt", "w") as out:
		out.write("Gene\tAssociated Disorder(s)\tSynopsis ID(s)\tNum Disorders\n")
		for gene in genesAndTheirDisorders:
			out.write("{}\t{}\t{}\t{}\n".format(gene, " | ".join(genesAndTheirDisorders[gene]["titles"]), " | ".join(genesAndTheirDisorders[gene]["indices"]), len(genesAndTheirDisorders[gene]["titles"])))
	return


# 2018-6-6
def getGeneListsForEachCategory():
	allSynopses = getAllSynopses()
	allCategoriesToGenes = defaultdict(set)

	print "Parsing synopses..."
	phenoMimsToGenesDict = makeMIMdict()
	for _, synopsis in allSynopses.iteritems():
		genes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict))
		for category in synopsis:
			if synopsis[category] == True:
				allCategoriesToGenes[category.replace("Exists", "")] |= genes

	for category in allCategoriesToGenes:
		with open(PROJECT + "geneListsForAllCategories/{date}/{cat}.txt".format(date = CURRENT_DATE_VERSION, cat = category), "w") as out:
			out.write("Gene lists for synopses that have a symptom in {cat}. Scraped {date}.\n".format(cat = category, date = CURRENT_DATE_VERSION))
			out.write("\n".join(allCategoriesToGenes[category]))



###################################################################################################

### Graphing bits; PLI.

###################################################################################################


def graphSetData():
	import pandas as pd
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt


	# plt.style.use('seaborn-deep')
	import seaborn as sns
	import sys

	sns.set_style('darkgrid')


	# disorderFileNames = ["iddData", "behavioralData", "kidneyData"]
	geneFileNames = ["iddGenes", "behavioralGenes", "uniqueIDDgenes", "uniqueBehavioralGenes", "iddAndBehavioralGenesIntersection", "kidneyGenes"]
	
	for geneList in geneFileNames:
		# data = pd.read_csv(PROJECT + geneList + ".txt", sep="\t")#, header=None)
		
		with open(PROJECT + geneList + ".txt", "r") as genes:
			header = genes.readline().strip().split("\t")
			inheritanceIndex = header.index("Inheritance")
			brainIndex = header.index("Brain Expressed?")
			inheritanceCounts = {}
			brainCounts = {}

			for gene in genes:
				gene = gene.strip().split("\t")
				inheritance = gene[inheritanceIndex]
				if inheritance not in inheritanceCounts:
					inheritanceCounts[inheritance] = 0
				inheritanceCounts[inheritance] += 1
				expressed = gene[brainIndex] 
				if expressed not in brainCounts:
					brainCounts[expressed] = 0
				brainCounts[expressed] += 1				



		inheritanceCounts = sorted(inheritanceCounts.iteritems(), key=lambda (k,v): (v,k), reverse=True)
		# A lot are one-offs; keep the real ones with more than one or two genes in them.
		inheritanceCounts = [(i, j) for (i,j) in inheritanceCounts if j > 10]
		inheritanceStrings = [i.replace(", ", ",\n") for (i,_) in inheritanceCounts]
		inheritanceValues = [j for (_,j) in inheritanceCounts]



		fig, (ax1, ax2) = plt.subplots(1, 2)

		# counts = data.Inheritance.value_counts()
		# inheritances = [i for i in counts.index]
		# inheritances = [",\n".join(inheritance.split(", ")) for inheritance in inheritances]

		# for name in names:
		# 	spaceIndex = name.rfind(" ")
		# 	if spaceIndex != -1:
		# 		spacedName = name[:spaceIndex] + "\n" + name[spaceIndex+1:]
		# 	names[names.index(name)] = spacedName
		# names = [name.replace(" ", "\n", 1) for name in names]
		# values = counts.values
		ax1.bar(range(len(inheritanceValues)), inheritanceValues)
		ax1.set_xticks(np.arange(len(inheritanceValues)))
		ax1.set_xticklabels(inheritanceStrings,rotation=70, fontsize=8)
		ax1.set_ylabel("Number of genes")
		ax1.set_title("Inheritance Method")

		# brain = data["Brain Expressed?"].value_counts()
		# brain.rename({"True": "Expressed", "False": "Not brain-expressed"})
		# names = brain.index
		# print "NAMES:", names
		# names = ["Expressed", "Not brain-expressed"]
		# values = brain.values	
		ax2.pie(brainCounts.values(), labels = brainCounts.keys(), autopct=lambda(x): '{:.0f}'.format(x * sum(brainCounts.values()) / 100.0))
		ax2.set_title("Brain Expression")
		ax2.axis('equal') # make the pie circular, not oval
		# fig.tight_layout(rect=[0, 0.03, 1, 0.95])
		plt.suptitle(geneList)
		plt.savefig(PROJECT + geneList + ".png", bbox_inches='tight')
		plt.close()


def graphPLIscores():

	pliScores = makeDictFromTwoColumns(PROJECT + "ReferenceFiles/pliScores.txt", "\t")


	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	# plt.style.use('seaborn-deep')
	import seaborn as sns
	import sys


	plt.figure(1)

	# # Original set comparison.
	# colors = ["red", "blue", "green"]
	# geneFileNames = ["iddGenes", "behavioralGenes", "kidneyGenes"]
	# labels = ["IDD Genes", "Behavioral Genes", "Kidney Genes"]
	# for geneFileName, color, label in zip(geneFileNames, colors, labels):

	# Graph pLI for all sets.
	FOLDER = PROJECT + "geneListsForAllCategories/"
	lowerLines = []
	upperLines = []
	for filename in os.listdir(FOLDER):

		geneList = makeSetFromFile(os.path.join(FOLDER, filename))
		name = filename.replace(".txt", "")
		# Format name nicely:
		setName = name[0].upper()
		for letter in name[1:]:
			if letter.istitle(): # check if capitalized
				setName += " "
			setName += letter
		notFound = []
		genesPLIscores = []
		for gene in geneList:
			try:
				genesPLIscores.append(float(pliScores[gene]))
			except KeyError as e:
				notFound.append(gene)
			except ValueError as e: # can't turn "NA" into float
				notFound.append(gene)
		# print "{} had {} (out of {} total) genes without pLI scores.".format(setName, len(notFound), len(geneList))
		
		# Empirical CDF: https://stackoverflow.com/a/11692365
		sortedScores = np.sort(genesPLIscores)
		yvals = np.arange(len(sortedScores))/float(len(sortedScores))
		# plt.plot(sortedScores, yvals, color=color, alpha=0.35, label=geneFileName.replace("Data",""))
		median =  np.median(sortedScores)

		if median < 0.0029:
			line, = plt.plot(sortedScores, yvals, alpha = 0.35, label = setName)
			lowerLines.append(line)
		if median > 0.33:
			line, = plt.plot(sortedScores, yvals, alpha = 0.35, label = setName)
			upperLines.append(line)
		# else:
		# 	plt.plot(sortedScores, yvals, alpha=0.35, label = None)

	# Create a legend for the first line.
	first_legend = plt.legend(handles=lowerLines, loc='lower right')
	# Add the legend manually to the current Axes.
	ax = plt.gca().add_artist(first_legend)

	# Create another legend for the second line.
	plt.legend(handles=upperLines, loc='upper left')

		# Original histogram
		# ax.hist(genesPLIscores, bins = 50, facecolor=color, alpha=0.35, label=geneFileName.replace("Data",""),\
			# ec='white', normed=True)
		# ax.legend(loc='upper right')

	plt.title("Gene Sets with Highest and Lowest Median PLI Scores")
	# plt.xscale('log')
	plt.xlabel('pLI Score')
	plt.ylabel('Percentage below this pLI value')
	# plt.legend(loc='lower right')
	# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.grid(True)

	plt.savefig(PROJECT + "highAndLowpLIscores.png", bbox_inches='tight')
	plt.close()



def findHowManyUniqueDisordersByPLI():

	genesToNumberDisorders = makeDictFromTwoColumns(PROJECT + "allGenes/genesToDisorders.txt", "\t")
	pliScores = makeDictFromTwoColumns(PROJECT + "ReferenceFiles/pliScores.txt", "\t")

	geneFileNames = ["iddGenesData", "behavioralGenesData", "kidneyGenesData"]
	for geneFileName in geneFileNames:
		print geneFileName
		geneDataList = makeSetFromFile(PROJECT + geneFileName + ".txt")
		geneList = [i.split("\t")[0] for i in geneDataList]
		
		notFound = []
		genesAndPLItuples = []
		for gene in geneList:
			try:
				genesAndPLItuples.append((gene,float(pliScores[gene])))
			except KeyError as e:
				notFound.append(gene)
			except ValueError as e: # can't turn "NA" into float
				notFound.append(gene)
		print "No pLI scores found for {} genes.".format(len(notFound))

		sortedScores = sorted(genesAndPLItuples, key=lambda (k,v): (v,k))
		topScores = [(gene,score) for (gene, score) in sortedScores if score > 0.9]

		totalGenes = len(topScores)
		sizeOfEachOf5Sets = totalGenes/5.0
		currentSetCount = 1
		numDisordersPerGene = []
		firstGene = None
		for i, (gene, _) in enumerate(topScores):
			if not firstGene:
				firstGene = gene
			try:
				numDisordersPerGene.append(float(genesToNumberDisorders[gene]))
				print "{}:{}".format(gene, float(genesToNumberDisorders[gene])),
			except KeyError as e:
				continue

			if i >= currentSetCount * sizeOfEachOf5Sets:
				print "Average disorders for set {} (n = {}; pLI = [{},{}]): {}".format(currentSetCount, len(numDisordersPerGene), pliScores[firstGene], pliScores[gene], float(sum(numDisordersPerGene))/len(numDisordersPerGene))
				numDisordersPerGene = []
				currentSetCount += 1
				firstGene = None
		print "\nAverage disorders for set {} (n = {}; pLI = [{},{}]): {}".format(currentSetCount, len(numDisordersPerGene), pliScores[firstGene], pliScores[gene], float(sum(numDisordersPerGene))/len(numDisordersPerGene))
		numDisordersPerGene = []
		currentSetCount += 1


		