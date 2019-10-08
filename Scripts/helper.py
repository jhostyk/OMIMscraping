# -*- coding: utf-8 -*-
### Helper functions and older versions of functions.

### Joseph Hostyk

### 2018-1-12

# PROJECT = "/nfs/labs/goldstein/jh3958/NeuroGenesPaper/"
PROJECT = "/Users/jhostyk/Dropbox/* Columbia/2a. Fall 2019/Symbolic Methods/Final Project/OMIMscraping/"
CURRENT_DATE_VERSION = "2019-08-02"
# import datetime
# CURRENT_DATE_VERSION = datetime.datetime.now().strftime("%Y-%m-%d")

from OMIMapiKey import *


import os
import sys
import time
import csv
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
	allSynopses = imp.load_source('', PROJECT + "Data/OMIMresults/{}_all.py".format(CURRENT_DATE_VERSION)).omim

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

# Lets us find associated genes more easily, later.
def makeMIMdict():

	phenoMimsToGenesDict = {}
	
	with open(PROJECT + "Data/ReferenceFiles/genemap2.txt", "r") as mims:
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

			for phenoMim in phenoMims:
				if phenoMim not in phenoMimsToGenesDict:
					phenoMimsToGenesDict[phenoMim] = []
				phenoMimsToGenesDict[phenoMim] += [(gene, geneMim)]

	return phenoMimsToGenesDict


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

def getOfficialGenes():

	officialGenes = set()

	hgncFile = PROJECT + "Data/ReferenceFiles/HGNCgenes.tsv"
	reader = csv.reader(open(hgncFile), delimiter='\t')
	header = next(reader)
	print header
	for entry in reader:
		entry = dict(zip(header, entry))
		officialGenes.add(entry["Approved symbol"])

	return officialGenes


def getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList):

	
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

	genes = [gene for gene in genes if gene in officialGeneList]

	return genes


def isKeyWordInAnySymptom(keyword, synopsis):

	for key in synopsis:
		if isinstance(synopsis[key], basestring) and keyword in synopsis[key].lower():
			return True
	if "oldFormat" in synopsis:
		return isKeyWordInAnySymptom(keyword, synopsis["oldFormat"])
	return False

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


# Don't know how to find genes' inheritance. Can only get per disorder.
# Want to see if for each gene, all the disorders they cause have the same type of inheritance.
def makeInheritanceFile():

	allSynopses = getAllSynopses()
	phenoMimsToGenesDict = makeMIMdict()
	brainGenes = makeSetFromFile(PROJECT + "ReferenceFiles/brainExpressedGenesFromDrSanders.txt")
	officialGeneList = getOfficialGenes()

	print "Parsing synopses..."
	genesInheritance = {}
	n = 0
	for _, synopsis in allSynopses.iteritems():
		sys.stdout.write("\rOn synopsis {}".format(n))
		sys.stdout.flush()
		n+=1
		genes, actualInheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict, officialGeneList)
		
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
	officialGeneList = getOfficialGenes()

	print "Parsing synopses..."
	genesData = {}
	dontWant = ["molecularBasis", "miscellaneous", "inheritance"]
	for _, synopsis in allSynopses.iteritems():

		genes, inheritance, _ = getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict, officialGeneList)
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


def getSynopsisInfo(synopsis, brainGenes, phenoMimsToGenesDict, officialGeneList):

	genes = getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList)

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
	
### Making reference files, using all the scraped data.
# Get and write each gene's associated disorders.
def getAllGenesDisorders():

	allSynopses = getAllSynopses()

	phenoMimsToGenesDict = makeMIMdict()
	officialGeneList = getOfficialGenes()
	

	genesAndTheirDisorders = {}

	print "Parsing synopses..."
	for synopsisID, synopsis in allSynopses.iteritems():

		genes = getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList)
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
	officialGeneList = getOfficialGenes()

	print "Parsing synopses..."
	phenoMimsToGenesDict = makeMIMdict()
	for _, synopsis in allSynopses.iteritems():
		genes = set(getGeneListFromSynopsis(synopsis, phenoMimsToGenesDict, officialGeneList))
		for category in synopsis:
			if synopsis[category] == True:
				allCategoriesToGenes[category.replace("Exists", "")] |= genes

	for category in allCategoriesToGenes:
		with open(PROJECT + "geneListsForAllCategories/{date}/{cat}.txt".format(date = CURRENT_DATE_VERSION, cat = category), "w") as out:
			out.write("Gene lists for synopses that have a symptom in {cat}. Scraped {date}.\n".format(cat = category, date = CURRENT_DATE_VERSION))
			out.write("\n".join(allCategoriesToGenes[category]))

