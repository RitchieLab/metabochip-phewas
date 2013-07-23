#!/usr/bin/env python

import collections
import math
import os
import sys


def go(dataFileName, snpIndexName, phenoIndexName, permu, snpListName, classListName, pvalCut):
	permuMax = 1024
	snpMax = 163840
	phenoMax = 416
	
	sys.stderr.write("reading SNP index from %s ..." % snpIndexName)
	snpN = {}
	snpName = {}
	with open(snpIndexName,'rU') as snpFile:
		line = snpFile.next().strip()
		if line != "snpNum,SNPID":
			raise Exception("unrecognized SNP index header: %s" % line)
		for line in snpFile:
			tok = line.strip().split(',',1)
			snpN[tok[1]] = int(tok[0])
			snpName[int(tok[0])] = tok[1]
		sys.stderr.write(" OK: %d SNPs\n" % len(snpN))
	
	sys.stderr.write("reading pheno index from %s ..." % phenoIndexName)
	phenoN = {}
	phenoName = {}
	phenoStudy = {}
	with open(phenoIndexName,'rU') as phenoFile:
		line = phenoFile.next().strip()
		if line != "phenoNum,phenotype,phenotype_transformation,study":
			raise Exception("unrecognized pheno index header: %s" % line)
		for line in phenoFile:
			tok = line.strip().split(',',1)
			phenoN[tok[1]] = int(tok[0])
			phenoName[int(tok[0])] = tok[1]
			phenoStudy[int(tok[0])] = tok[1].split(',')[-1]
		sys.stderr.write(" OK: %d pheno,xform,study tuples\n" % len(phenoN))
	
	sys.stderr.write("reading class list from %s ..." % classListName)
	phenoNameClass = {}
	with open(classListName,'rU') as classFile:
		line = classFile.next().strip()
		if line != "phenotype,class":
			raise Exception("unrecognized class list header: %s" % line)
		for line in classFile:
			tok = line.strip().split(',',1)
			if tok[1] != "" and tok[1] != "NULL":
				phenoNameClass[tok[0].lower()] = tok[1]
		sys.stderr.write(" OK: %d phenotype class assignments\n" % len(phenoNameClass))
	
	phenoClass = {}
	phenoNoClass = set()
	for phenotag in phenoN:
		phenotype = phenotag.split(',',1)[0].lower()
		if phenotype in phenoNameClass:
			phenoClass[phenoN[phenotag]] = phenoNameClass[phenotype]
		else:
			phenoClass[phenoN[phenotag]] = None
			phenoNoClass.add(phenotype)
	if len(phenoNoClass) > 0:
		sys.stderr.write("WARNING: %d phenotypes missing from class list:\n  %s\n" % (len(phenoNoClass), ",".join(sorted(phenoNoClass))))
	
	if permu < 0 or permu >= permuMax:
		raise Exception("permutation %d out of bounds (0-%d)\n" % (permu,permuMax))
	
	sys.stderr.write("reading SNP list from %s ..." % snpListName)
	snpList = []
	snpSkip = []
	with open(snpListName,'rU') as snpFile:
		line = snpFile.next().strip()
		if line != "SNPID":
			raise Exception("unrecognized SNP list header: %s" % line)
		for line in snpFile:
			snp = line.strip()
			if snp in snpN:
				snpList.append(snpN[snp])
			else:
				snpSkip.append(snp)
	snpList.sort()
	sys.stderr.write(" OK: %d SNPs\n" % len(snpList))
	if len(snpSkip) > 0:
		sys.stderr.write("WARNING: %d SNPs not recognized:\n  %s\n" % (len(snpSkip), ",".join(snpSkip)))
	
	phenoList = []
	for p in phenoName:
		if p in phenoClass and phenoClass[p]:
			phenoList.append(p)
	phenoList.sort()
	
	sys.stderr.write(
"""analyzing %s ...
  1 permutation: %d
  %d SNP%s
  %d phenotype%s
""" % (
		dataFileName,
		permu,
		len(snpList), (": %s" % snpName[snpList[0]]) if len(snpList) == 1 else "s",
		len(phenoList), (": %s" % phenoName[phenoList[0]]) if len(phenoList) == 1 else "s",
	))
	
	numMissing = numMin = numMax = numNA = numLow = numHigh = 0
	snpsR1 = snpsR2 = snpsR3 = 0
	with open(dataFileName,'rb') as dataFile:
		for s in snpList:
			classStudies = collections.defaultdict(set)
			pos = ((((permu * snpMax) + s) * phenoMax) + 0) * 3
			dataFile.seek(pos)
			pdata = dataFile.read(phenoMax * 3)
			for p in phenoList:
				data = pdata[3*p : 3*p+3]
				if data == "" or data == "\xFF\xFF\xFF": #never written
					pval = None
					numMissing += 1
				elif data == "\x00\x00\xFF": #NA
					pval = None
					numNA += 1
				elif data == "\x00\x01\xFF": #0.001e-255 lower limit
					pval = 0.001e-255
					numMin += 1
				elif data == "\xFF\xFF\x00": #65.535e-0 upper limit
					pval = 65.535
					numMax += 1
				else:
					pval = (ord(data[0])*256.0 + ord(data[1])) / 1000.0 * math.pow(10, -ord(data[2]))
					if pval <= pvalCut:
						numLow += 1
					else:
						numHigh += 1
				if pval != None:
					if pval <= pvalCut:
						classStudies[phenoClass[p]].add(phenoStudy[p])
			#foreach pheno
			n = sum((1 if len(classStudies[c]) > 1 else 0) for c in classStudies)
			if n >= 1:
				snpsR1 += 1
				if n >= 2:
					snpsR2 += 1
					if n >= 3:
						snpsR3 += 1
		#foreach snp
	#with dataFile
	sys.stderr.write("... done.\n")
	
	numTotal = numMissing+numMin+numMax+numNA+numLow+numHigh
	
	sys.stderr.write(
"""out of %d possible data points:
  %d (%1.2f%%) missing from source data
  %d (%1.2f%%) marked 'NA' in source data
  %d (%1.2f%%) above %g, including %d (%1.2f%%) out-of-range
  %d (%1.2f%%) below %g, including %d (%1.2f%%) out-of-range
out of %d SNPs:
  %d (%1.2f%%) replicate at least 1 class across 2 or more studies
  %d (%1.2f%%) replicate at least 2 classes across 2 or more studies
  %d (%1.2f%%) replicate at least 3 classes across 2 or more studies
""" % (
		numTotal,
		numMissing, 100.0*numMissing/numTotal,
		numNA, 100.0*numNA/numTotal,
		numMax+numHigh, 100.0*(numMax+numHigh)/numTotal, pvalCut, numMax, 100.0*numMax/numTotal,
		numMin+numLow, 100.0*(numMin+numLow)/numTotal, pvalCut, numMin, 100.0*numMin/numTotal,
		len(snpList),
		snpsR1, 100.0*snpsR1/len(snpList),
		snpsR2, 100.0*snpsR2/len(snpList),
		snpsR3, 100.0*snpsR3/len(snpList)
	))
	sys.stdout.write("%d\t%d\t%d\t%d\n" % (permu,snpsR1,snpsR2,snpsR3))
#go()


if __name__ == "__main__":
	if len(sys.argv) < 8:
		sys.stderr.write(
"""Usage: %s <datafile> <snpindex> <phenoindex> <permu> <snplist> <classlist> <pvalcutoff>

The <datafile> will be scanned using the <snpindex> and <phenoindex> to
generate phenotype class replication statistics for the given <permu> based on
the given <pvalcutoff>.  Only SNPs appearing in the <snplist> file and
phenotypes with classes in the <classlist> file will be considered.
""" % os.path.basename(sys.argv[0]))
		exit()
	
	dataFileName = sys.argv[1]
	snpIndexName = sys.argv[2]
	phenoIndexName = sys.argv[3]
	permu = int(sys.argv[4])
	snpListName = sys.argv[5]
	classListName = sys.argv[6]
	pvalCut = float(sys.argv[7])
	
	go(dataFileName, snpIndexName, phenoIndexName, permu, snpListName, classListName, pvalCut)
#__main__
