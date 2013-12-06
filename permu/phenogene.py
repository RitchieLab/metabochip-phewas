#!/usr/bin/env python

"""
/* for each study/pheno/xform/gene, identify the best SNP association within 10k of the gene */
SELECT
  d.study, d.substudy, d.phenotype, d.phenotype_transformation, ni.entrezID,
  SUBSTRING_INDEX(GROUP_CONCAT(d.id ORDER BY d.pval),',',1) AS id
FROM ncbi_geneinfo AS ni
JOIN ncbi_gene2refseq AS nr USING (entrezID)
JOIN mcpwas_markers AS m
  ON m.chrom = ni.chr AND m.pos BETWEEN nr.startPos-10000 AND nr.endPos+10000
JOIN mcpwas_data AS d
  ON d.SNPID = m.id AND d.pval > 0 AND d.AF_coded_all >= 0.01
GROUP BY
  d.study, d.substudy, d.phenotype, d.phenotype_transformation, ni.entrezID
"""

import sys
import MySQLdb
import Image

study = "MEC"
geneCutoff = 0
phenoCutoff = 0
if len(sys.argv) >= 2:
	study = sys.argv[1]
	if len(sys.argv) >= 3:
		geneCutoff = float(sys.argv[2])
		if len(sys.argv) >= 4:
			phenoCutoff = float(sys.argv[3])
sys.stderr.write("""tabulating gene-phenotype/transform data
	study: %s
	gene cutoff: %f
	pheno cutoff: %f
""" % (study,geneCutoff,phenoCutoff))
	
try:
	try:
		# connect
		conn = MySQLdb.connect(
			host = 'badger',
			user = 'atf3',
			passwd = 'patO9FTPXUV0JonedaLe1COO',
			db = 'ritchie_lab'
		)
		curs = conn.cursor()
		
		# load relevant genes
		numGene = 0
		geneTotal = {}
		geneBest = {}
		geneName = {}
		curs.execute("SELECT DISTINCT entrezID,symbol FROM mcpwas_phenogene JOIN ncbi_geneinfo USING (entrezID) WHERE study = %s", (study))
		row = curs.fetchone()
		while (row):
			numGene += 1
			geneTotal[row[0]] = 0
			geneBest[row[0]] = 0
			geneName[row[0]] = row[1]
			row = curs.fetchone()
		sys.stderr.write("loaded %d genes\n" % (numGene))
		
		# load relevant phenotypes
		numPheno = 0
		numPhenoXform = 0
		phenoXformTotal = {}
		phenoXformBest = {}
		curs.execute("SELECT DISTINCT LOWER(phenotype),phenotype_transformation FROM mcpwas_phenogene WHERE study = %s", (study))
		row = curs.fetchone()
		while (row):
			if (row[0] not in phenoXformTotal):
				numPheno += 1
				phenoXformTotal[row[0]] = {}
				phenoXformBest[row[0]] = {}
			numPhenoXform += 1
			phenoXformTotal[row[0]][row[1]] = 0
			phenoXformBest[row[0]][row[1]] = 0
			row = curs.fetchone()
		sys.stderr.write("loaded %d phenotype-transforms (%d phenotypes)\n" % (numPhenoXform, numPheno))
		
		# load gene-pheno-xform data
		numAssoc = 0
		genePhenoXformVal = {}
		curs.execute("SELECT pg.entrezID,LOWER(pg.phenotype),pg.phenotype_transformation,-log10(d.pval) FROM mcpwas_phenogene AS pg JOIN mcpwas_data AS d USING (id) WHERE pg.study = %s", (study))
		row = curs.fetchone()
		minVal = maxVal = row[3]
		while (row):
			numAssoc += 1
			if (row[0] not in genePhenoXformVal):
				genePhenoXformVal[row[0]] = {row[1]: {}}
			elif (row[1] not in genePhenoXformVal[row[0]]):
				genePhenoXformVal[row[0]][row[1]] = {}
			genePhenoXformVal[row[0]][row[1]][row[2]] = row[3]
			geneTotal[row[0]] += row[3]
			geneBest[row[0]] = max(geneBest[row[0]], row[3])
			phenoXformTotal[row[1]][row[2]] += row[3]
			phenoXformBest[row[1]][row[2]] = max(phenoXformBest[row[1]][row[2]], row[3])
			minVal = min(minVal, row[3])
			maxVal = max(maxVal, row[3])
			row = curs.fetchone()
		sys.stderr.write("loaded %d gene-phenotype-transform tuples out of %d possible (-log10(pval) = %f ... %f)\n" % (numAssoc, numGene*numPhenoXform, minVal, maxVal))
		
		# sort phenotypes and genes by total score
		phenoTotal = {}
		for pheno in phenoXformTotal:
			phenoTotal[pheno] = max(phenoXformTotal[pheno].values())
		phenoOrder = sorted(phenoTotal, key=phenoTotal.get, reverse=True)
		geneOrder = sorted(geneTotal, key=geneTotal.get, reverse=True)
		sys.stderr.write("best phenotypes:\n\t%s (%f)\n\t%s (%f)\n\t%s (%f)\n" % (
				phenoOrder[0], phenoTotal[phenoOrder[0]] / numGene,
				phenoOrder[1], phenoTotal[phenoOrder[1]] / numGene,
				phenoOrder[2], phenoTotal[phenoOrder[2]] / numGene
		))
		sys.stderr.write("best genes:\n\t%s (%f)\n\t%s (%f)\n\t%s (%f)\n" % (
				geneName[geneOrder[0]], geneTotal[geneOrder[0]] / numPhenoXform,
				geneName[geneOrder[1]], geneTotal[geneOrder[1]] / numPhenoXform,
				geneName[geneOrder[2]], geneTotal[geneOrder[2]] / numPhenoXform
		))
		
		# output associations as CSV
		sys.stdout.write('gene')
		skipGene = 0
		skipPheno = 0
		for pheno in phenoOrder:
			for xform in phenoXformBest[pheno]:
				if phenoXformBest[pheno][xform] < phenoCutoff:
					skipPheno += 1
				else:
					sys.stdout.write(',"%s (%s)"' % (pheno,xform))
		sys.stdout.write("\n")
		for gene in geneOrder:
			if geneBest[gene] < geneCutoff:
				skipGene += 1
			else:
				sys.stdout.write('"%s:%d"' % (geneName[gene],gene))
				for pheno in phenoOrder:
					for xform in phenoXformBest[pheno]:
						if phenoXformBest[pheno][xform] < phenoCutoff:
							pass
						elif (gene in genePhenoXformVal) and (pheno in genePhenoXformVal[gene]) and (xform in genePhenoXformVal[gene][pheno]):
							sys.stdout.write(",%f" % genePhenoXformVal[gene][pheno][xform])
						else:
							sys.stdout.write(",NA")
				sys.stdout.write("\n")
		sys.stderr.write("wrote %dx%d gene/pheno table (dropped %d gene rows, %d pheno-xform cols)\n" % (
				numGene-skipGene, numPhenoXform-skipPheno, skipGene, skipPheno
		))
		
	except MySQLdb.Error, e:
		print "MySQL Error #%d: %s" % (e.args[0], e.args[1])
		sys.exit(1)
	
finally:
	if (conn):
		conn.close()
