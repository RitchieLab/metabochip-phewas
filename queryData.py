#!/usr/bin/env python

import sys
import os
import math


def go(dataFileName, snpIndexName, phenoIndexName, pvalCut, permuList, snpList, phenoList):
	permuMax = 1024
	snpMax = 163840
	phenoMax = 416
	
	sys.stdout.write("reading SNP index from %s ..." % snpIndexName)
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
		sys.stdout.write(" OK: %d SNPs\n" % len(snpN))
	
	sys.stdout.write("reading pheno index from %s ..." % phenoIndexName)
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
		sys.stdout.write(" OK: %d pheno,xform,study tuples\n" % len(phenoN))
	
	permuUsed = 1000
	snpUsed = len(snpN) - 1
	phenoUsed = len(phenoN) - 1
	
	# validate permutation list
	if permuList:
		permuListN = []
		for permu in permuList:
			p = int(permu)
			if p < 0 or p >= permuMax:
				sys.stderr.write("permutation %d out of bounds (0-%d)\n" % (p,permuMax))
				sys.exit(1)
			permuListN.append(p)
		permuListN.sort()
		permuList = permuListN
	else:
		permuList = range(0,permuMax)
	
	# validate SNP list
	if snpList:
		snpListN = []
		for snp in snpList:
			if snp not in snpN:
				sys.stderr.write("unrecognized SNP: %s\n" % snp)
				sys.exit(1)
			snpListN.append(snpN[snp])
		snpListN.sort()
		snpList = snpListN
	else:
		snpList = range(0,len(snpN))
	
	# validate pheno list
	if phenoList:
		phenoListN = []
		for pheno in phenoList:
			if pheno not in phenoN:
				sys.stderr.write("unrecognized phenotype: %s\n" % pheno)
				sys.exit(1)
			phenoListN.append(phenoN[pheno])
		phenoListN.sort()
		phenoList = phenoListN
	else:
		phenoList = range(0,len(phenoN))
	
	sys.stdout.write(
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
	numPad = numMissing = numMin = numMax = numNA = numLow = numHigh = 0
	
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
					elif data == "\x00\x01\xFF": #0.001e-255 limit
						numMin += 1
					elif data == "\xFF\xFF\x00": #65.535e-0 out-of-range
						numMax += 1
					elif data == "\x00\x00\xFF": #NA
						numNA += 1
					#elif (ord(data[0])*256.0 + ord(data[1])) * math.pow(10, -ord(data[2])) <= pvalCut:
					elif (ord(data[0])*256.0 + ord(data[1])) / 1000.0 * math.pow(10, -ord(data[2])) <= pvalCut:
						numLow += 1
					else:
						numHigh += 1
				#foreach pheno
			#foreach snp
		#foreach permu
	#with dataFile
	sys.stdout.write("... done.\n")
	
	numTotal = numPad+numMissing+numMin+numMax+numNA+numLow+numHigh
	sys.stdout.write(
"""out of %d possible data points:
  %d (%1.2f%%) in unused padding regions
  %d (%1.2f%%) missing from source data
  %d (%1.2f%%) marked 'NA' in source data
  %d (%1.2f%%) above %g, including %d (%1.2f%%) out-of-range
  %d (%1.2f%%) below %g, including %d (%1.2f%%) out-of-range
""" % (
		numTotal,
		numPad, 100.0*numPad/numTotal,
		numMissing, 100.0*numMissing/numTotal,
		numNA, 100.0*numNA/numTotal,
		numMax+numHigh, 100.0*(numMax+numHigh)/numTotal, pvalCut, numMax, 100.0*numMax/numTotal,
		numMin+numLow, 100.0*(numMin+numLow)/numTotal, pvalCut, numMin, 100.0*numMin/numTotal
	))
#go()


if __name__ == "__main__":
	if len(sys.argv) < 5:
		sys.stderr.write(
"""Usage: %s <datafile> <snpindex> <phenoindex> <pvalCutoff> [permu] [snp] [pheno]

The <datafile> will be scanned using the <snpindex> and <phenoindex> to
generate statistics based on the given <pvalCutoff>.

The search may optionally be limited to one or more [permu]tations, [snp]s or
[pheno]types.  Multiple values may be joined with a "+" (i.e.
rs974059+rs1181871), and all values may be specified with a "-" (in order to
specify a later argument).
""" % os.path.basename(sys.argv[0]))
		exit()
	
	dataFileName = sys.argv[1]
	snpIndexName = sys.argv[2]
	phenoIndexName = sys.argv[3]
	pvalCut = float(sys.argv[4])
	permuList = snpList = phenoList = None
	
	if len(sys.argv) >= 6 and sys.argv[5] != "-":
		permuList = sys.argv[5].split('+')
	if len(sys.argv) >= 7 and sys.argv[6] != "-":
		snpList = sys.argv[6].split('+')
	if len(sys.argv) >= 8 and sys.argv[7] != "-":
		phenoList = sys.argv[7].split('+')
	
	go(dataFileName, snpIndexName, phenoIndexName, pvalCut, permuList, snpList, phenoList)
#__main__
