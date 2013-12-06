#!/usr/bin/env python

import sys
import os
import math


def go(dataFileName, snpIndexName, phenoIndexName, phenoList):
	permuMax = 1024
	snpMax = 163840
	phenoMax = 416
	
	sys.stderr.write("reading SNP index from %s ..." % snpIndexName)
	snpN = {}
	snpName = {}
	with open(snpIndexName,'r') as snpFile:
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
	with open(phenoIndexName,'r') as phenoFile:
		line = phenoFile.next().strip()
		if line != "phenoNum,phenotype,phenotype_transformation,study":
			raise Exception("unrecognized pheno index header: %s" % line)
		for line in phenoFile:
			tok = line.strip().split(',',1)
			phenoN[tok[1]] = int(tok[0])
			phenoName[int(tok[0])] = tok[1]
	sys.stderr.write(" OK: %d pheno,xform,study tuples\n" % len(phenoN))
	
	permuUsed = 1000
	snpUsed = len(snpN) - 1
	phenoUsed = len(phenoN) - 1
	
	permuList = range(0,permuUsed)
	snpList = range(0,snpUsed+1)
	
	# validate pheno list
	phenoListN = []
	for pheno in phenoList:
		if pheno not in phenoN:
			sys.stderr.write("unrecognized phenotype: %s\n" % pheno)
			sys.exit(1)
		phenoListN.append(phenoN[pheno])
	phenoListN.sort()
	phenoList = phenoListN
	
	sys.stderr.write(
"""analyzing %s ...
  %d permutation%s
  %d SNP%s
  %d phenotype%s
""" % (
		dataFileName,
		len(permuList), (": %d" % permuList[0]) if len(permuList) == 1 else "s",
		len(snpList), (": %s" % snpName[snpList[0]]) if len(snpList) == 1 else "s",
		len(phenoList), (": %s" % phenoName[phenoList[0]]) if len(phenoList) == 1 else "s",
	))
	
	numPad = numMissing = numNA = numMax = numMin = numVal = 0
	snpPvalMax = { s:-9 for s in snpList }
	snpPvalMin = { s:9 for s in snpList }
	snpPvalSum = { s:0 for s in snpList }
	snpPvalNum = { s:0 for s in snpList }
	with open(dataFileName,'rb') as dataFile:
		for r in permuList:
			if r > permuUsed:
				numPad += len(snpList) * len(phenoList)
				continue
			for s in snpList:
				if s > snpUsed:
					numPad += len(phenoList)
					continue
				for p in phenoList:
					if p > phenoUsed:
						numPad += 1
						continue
					pos = ((((r * snpMax) + s) * phenoMax) + p) * 3
					dataFile.seek(pos)
					data = dataFile.read(3)
					if data == "" or data == "\xFF\xFF\xFF": #never written
						numMissing += 1
					elif data == "\x00\x00\xFF": #NA
						numNA += 1
					elif data == "\xFF\xFF\x00": #65.535e-0 out-of-range
						numMax += 1
						snpPvalMax[s] = 65.535
						snpPvalSum[s] += 65.535
						snpPvalNum[s] += 1
					elif data == "\x00\x01\xFF": #0.001e-255 limit
						numMin += 1
						snpPvalMin[s] = 0.001e-255
						snpPvalSum[s] += 0.001e-255
						snpPvalNum[s] += 1
					else:
						numVal += 1
						pvalue = (ord(data[0])*256.0 + ord(data[1])) / 1000.0 * math.pow(10, -ord(data[2]))
						snpPvalMin[s] = min(snpPvalMin[s], pvalue)
						snpPvalMax[s] = max(snpPvalMax[s], pvalue)
						snpPvalSum[s] += pvalue
						snpPvalNum[s] += 1
				#foreach pheno
			#foreach snp
		#foreach permu
	#with dataFile
	sys.stderr.write("... done.\n")
	
	numTotal = numPad+numMissing+numNA+numMax+numMin+numVal
	sys.stderr.write(
"""out of %d possible data points:
  %d (%1.2f%%) in unused padding regions
  %d (%1.2f%%) missing from source data
  %d (%1.2f%%) marked 'NA' in source data
  %d (%1.2f%%) out-of-range high (65.535)
  %d (%1.2f%%) out-of-range low (1e-258)
  %d (%1.2f%%) tallied
""" % (
		numTotal,
		numPad, 100.0*numPad/numTotal,
		numMissing, 100.0*numMissing/numTotal,
		numNA, 100.0*numNA/numTotal,
		numMax, 100.0*numMax/numTotal,
		numMin, 100.0*numMin/numTotal,
		numVal, 100.0*numVal/numTotal
	))
	
	sys.stdout.write("phenotype\tphenotype_transformation\tstudy\tid\tnum\tmin\tavg\tmax\n")
	for s in snpList:
		if snpPvalNum[s] <= 0:
			sys.stdout.write("%s\t%s\t%d\tNULL\tNULL\tNULL\n" % (
					phenoName[phenoList[0]].replace(",","\t"),
					snpName[s],
					snpPvalNum[s]
			))
		else:
			sys.stdout.write("%s\t%s\t%d\t%g\t%g\t%g\n" % (
					phenoName[phenoList[0]].replace(",","\t"),
					snpName[s],
					snpPvalNum[s],
					snpPvalMin[s],
					snpPvalSum[s]/snpPvalNum[s],
					snpPvalMax[s]
			))
#go()


if __name__ == "__main__":
	if len(sys.argv) < 5:
		sys.stderr.write(
"""Usage: %s <datafile> <snpindex> <phenoindex> <pheno>

The <datafile> will be scanned using the <snpindex> and <phenoindex> to
generate statistics for the given <pheno>.
""" % os.path.basename(sys.argv[0]))
		exit()
	
	dataFileName = sys.argv[1]
	snpIndexName = sys.argv[2]
	phenoIndexName = sys.argv[3]
	phenoList = sys.argv[4].split('+')
	
	go(dataFileName, snpIndexName, phenoIndexName, phenoList)
#__main__
