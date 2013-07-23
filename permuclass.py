#!/usr/bin/env python

import sys
import os
import math


def go(dataFileName, snpIndexName, phenoIndexName, classListName, pvalCut, permu, snpListName, study):
	permuMax = 1024
	snpMax = 163840
	phenoMax = 416
	
	sys.stdout.write("reading SNP index from %s ..." % snpIndexName)
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
		sys.stdout.write(" OK: %d SNPs\n" % len(snpN))
	
	sys.stdout.write("reading pheno index from %s ..." % phenoIndexName)
	phenoN = {}
	phenoName = {}
	with open(phenoIndexName,'rU') as phenoFile:
		line = phenoFile.next().strip()
		if line != "phenoNum,phenotype,phenotype_transformation,study":
			raise Exception("unrecognized pheno index header: %s" % line)
		for line in phenoFile:
			tok = line.strip().split(',',1)
			phenoN[tok[1]] = int(tok[0])
			phenoName[int(tok[0])] = tok[1]
		sys.stdout.write(" OK: %d pheno,xform,study tuples\n" % len(phenoN))
	
	sys.stdout.write("reading class list from %s ..." % classListName)
	phenoNameClass = {}
	with open(classListName,'rU') as classFile:
		line = classFile.next().strip()
		if line != "phenotype,class":
			raise Exception("unrecognized class list header: %s" % line)
		for line in classFile:
			tok = line.strip().split(',',1)
			if tok[1] != "" and tok[1] != "NULL":
				phenoNameClass[tok[0].lower()] = tok[1]
		sys.stdout.write(" OK: %d phenotype class assignments\n" % len(phenoNameClass))
	
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
		sys.stderr.write("permutation %d out of bounds (0-%d)\n" % (permu,permuMax))
		sys.exit(1)
	
	sys.stdout.write("reading SNP list from %s ..." % snpListName)
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
	sys.stdout.write(" OK: %d SNPs\n" % len(snpList))
	if len(snpSkip) > 0:
		sys.stderr.write("WARNING: %d SNPs not recognized:\n  %s\n" % (len(snpSkip), ",".join(snpSkip)))
	
	phenoList = []
	for p in phenoName:
		if p in phenoClass and phenoClass[p]:
			if (not study) or phenoName[p].endswith(",%s" % study):
				phenoList.append(p)
	phenoList.sort()
	
	sys.stdout.write(
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
	
	numSkip = numMissing = numMin = numMax = numNA = numLow = numHigh = 0
	bestPval = 9.9
	bestS = bestP = None
	snpC1 = snpC2 = snpC3 = snpC4 = snpC5 = snpC6 = 0
	with open(dataFileName,'rb') as dataFile:
		for s in snpList:
			if pvalCut >= 0:
				classHit = set()
			pos = ((((permu * snpMax) + s) * phenoMax) + 0) * 3
			dataFile.seek(pos)
			pdata = dataFile.read(phenoMax * 3)
			for p in phenoList:
				if p in phenoClass and phenoClass[p]:
					data = pdata[3*p : 3*p+3]
					if data == "" or data == "\xFF\xFF\xFF": #never written
						numMissing += 1
					elif data == "\x00\x01\xFF": #0.001e-255 limit
						if pvalCut < 0 and numMin == 0:
							bestPval = 0.001e-255
							bestS = s
							bestP = p
						numMin += 1
					elif data == "\xFF\xFF\x00": #65.535e-0 out-of-range
						numMax += 1
					elif data == "\x00\x00\xFF": #NA
						numNA += 1
					else:
						pval = (ord(data[0])*256.0 + ord(data[1])) / 1000.0 * math.pow(10, -ord(data[2]))
						if pvalCut < 0:
							if pval < bestPval and numMin == 0:
								bestPval = pval
								bestS = s
								bestP = p
							numLow += 1
						elif pval <= pvalCut:
							classHit.add(phenoClass[p])
							numLow += 1
						else:
							numHigh += 1
				else:
					numSkip += 1
			#foreach pheno
			if pvalCut >= 0:
				if len(classHit) >= 1:
					snpC1 += 1
					if len(classHit) >= 2:
						snpC2 += 1
						if len(classHit) >= 3:
							snpC3 += 1
							if len(classHit) >= 4:
								snpC4 += 1
								if len(classHit) >= 5:
									snpC5 += 1
									if len(classHit) >= 6:
										snpC6 += 1
		#foreach snp
	#with dataFile
	sys.stdout.write("... done.\n")
	
	numTotal = numSkip+numMissing+numMin+numMax+numNA+numLow+numHigh
	
	if pvalCut < 0:
		sys.stdout.write(
"""out of %d possible data points:
  %d (%1.2f%%) skipped due to class or study filter
  %d (%1.2f%%) missing from source data
  %d (%1.2f%%) marked 'NA' in source data
  %d below .001e-255 limit
  best pval/SNP/pheno:
%1.6e
%s
%s
""" % (
			numTotal,
			numSkip, 100.0*numSkip/numTotal,
			numMissing, 100.0*numMissing/numTotal,
			numNA, 100.0*numNA/numTotal,
			numMin,
			bestPval,
			snpName[bestS],
			phenoName[bestP]
		))
	else:
		sys.stdout.write(
"""out of %d possible data points:
  %d (%1.2f%%) skipped due to missing phenotype class
  %d (%1.2f%%) missing from source data
  %d (%1.2f%%) marked 'NA' in source data
  %d (%1.2f%%) above %g, including %d (%1.2f%%) out-of-range
  %d (%1.2f%%) below %g, including %d (%1.2f%%) out-of-range
out of %d SNPs:
  %d (%1.2f%%) associated with at least 1 phenotype class
  %d (%1.2f%%) associated with at least 2 different classes
  %d (%1.2f%%) associated with at least 3 different classes
  %d (%1.2f%%) associated with at least 4 different classes
  %d (%1.2f%%) associated with at least 5 different classes
  %d (%1.2f%%) associated with at least 6 different classes
""" % (
			numTotal,
			numSkip, 100.0*numSkip/numTotal,
			numMissing, 100.0*numMissing/numTotal,
			numNA, 100.0*numNA/numTotal,
			numMax+numHigh, 100.0*(numMax+numHigh)/numTotal, pvalCut, numMax, 100.0*numMax/numTotal,
			numMin+numLow, 100.0*(numMin+numLow)/numTotal, pvalCut, numMin, 100.0*numMin/numTotal,
			len(snpList),
			snpC1, 100.0*snpC1/len(snpList),
			snpC2, 100.0*snpC2/len(snpList),
			snpC3, 100.0*snpC3/len(snpList),
			snpC4, 100.0*snpC4/len(snpList),
			snpC5, 100.0*snpC5/len(snpList),
			snpC6, 100.0*snpC6/len(snpList)
		))
#go()


if __name__ == "__main__":
	if len(sys.argv) < 8:
		sys.stderr.write(
"""Usage: %s <datafile> <snpindex> <phenoindex> <classlist> <pvalcutoff> <permu> <snplist> [study]

The <datafile> will be scanned using the <snpindex> and <phenoindex> to
generate phenotype class association statistics for the given <permu> based on
the given <pvalcutoff>.  Only SNPs appearing in the <snplist> file will be
considered.
""" % os.path.basename(sys.argv[0]))
		exit()
	
	dataFileName = sys.argv[1]
	snpIndexName = sys.argv[2]
	phenoIndexName = sys.argv[3]
	classListName = sys.argv[4]
	pvalCut = float(sys.argv[5])
	permu = int(sys.argv[6])
	snpListName = sys.argv[7]
	study = sys.argv[8] if len(sys.argv) > 8 else None
	
	go(dataFileName, snpIndexName, phenoIndexName, classListName, pvalCut, permu, snpListName, study)
#__main__
