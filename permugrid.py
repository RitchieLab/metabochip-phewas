#!/usr/bin/env python

import sys
import os
import math


def go(dataFileName, snpIndexName, phenoIndexName, classListName, snpListName):
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
	with open(phenoIndexName,'rU') as phenoFile:
		line = phenoFile.next().strip()
		if line != "phenoNum,phenotype,phenotype_transformation,study":
			raise Exception("unrecognized pheno index header: %s" % line)
		for line in phenoFile:
			tok = line.strip().split(',',1)
			phenoN[tok[1]] = int(tok[0])
			phenoName[int(tok[0])] = tok[1]
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
	
	permuList = range(0,permuMax+1)
	
	phenoList = []
	for p in phenoClass:
		if phenoClass[p]:
			phenoList.append(p)
	phenoList.sort()
	
	sys.stderr.write(
"""analyzing %s ...
  all permutations
  %d SNP%s
  %d phenotype%s
""" % (
		dataFileName,
		len(snpList), (": %s" % snpName[snpList[0]]) if len(snpList) == 1 else "s",
		len(phenoList), (": %s" % phenoName[phenoList[0]]) if len(phenoList) == 1 else "s",
	))
	
	best = {}
	for s in snpList:
		best[s] = {}
		for p in phenoList:
			best[s][p] = 9
	with open(dataFileName,'rb') as dataFile:
		for r in permuList:
			for s in snpList:
				pos = ((((r * snpMax) + s) * phenoMax) + 0) * 3
				dataFile.seek(pos)
				pdata = dataFile.read(phenoMax * 3)
				for p in phenoList:
					data = pdata[3*p : 3*p+3]
					if data == "" or data == "\xFF\xFF\xFF": #never written
						pass
					elif data == "\x00\x01\xFF": #0.001e-255 limit
						best[s][p] = 0.001e-255
					elif data == "\xFF\xFF\x00": #65.535e-0 out-of-range
						pass
					elif data == "\x00\x00\xFF": #NA
						pass
					else:
						pval = (ord(data[0])*256.0 + ord(data[1])) / 1000.0 * math.pow(10, -ord(data[2]))
						best[s][p] = min(best[s][p], pval)
				#foreach pheno
			#foreach snp
		#foreach permu
	#with dataFile
	sys.stderr.write("... done.\n")
	
	sys.stdout.write("SNPID")
	for p in phenoList:
		sys.stdout.write(",\"%s\"" % phenoName[p])
	sys.stdout.write("\n")
	for s in snpList:
		sys.stdout.write(snpName[s])
		for p in phenoList:
			sys.stdout.write(",%1.3e" % best[s][p])
		sys.stdout.write("\n")
#go()


if __name__ == "__main__":
	if len(sys.argv) < 6:
		sys.stderr.write(
"""Usage: %s <datafile> <snpindex> <phenoindex> <classlist> <snplist>

The <datafile> will be scanned using the <snpindex> and <phenoindex> to
generate a grid of best p-values.  Only phenotypes with a class in <classlist>
and SNPs appearing in the <snplist> file will be considered.
""" % os.path.basename(sys.argv[0]))
		exit()
	
	dataFileName = sys.argv[1]
	snpIndexName = sys.argv[2]
	phenoIndexName = sys.argv[3]
	classListName = sys.argv[4]
	snpListName = sys.argv[5]
	
	go(dataFileName, snpIndexName, phenoIndexName, classListName, snpListName)
#__main__
