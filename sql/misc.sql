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


/* count NCBI genes which overlap at least one other gene */
SELECT STRAIGHT_JOIN
  COUNT(DISTINCT r1.entrezID)
FROM ncbi_gene2refseq AS r1
JOIN ncbi_gene2refseq AS r2
  ON r2.geneAcc LIKE 'NC%' AND r2.startPos >= r1.startPos AND r2.startPos <= r1.endPos AND r2.entrezID != r1.entrezID
JOIN ncbi_geneinfo AS i1
  ON i1.entrezID = r1.entrezID
JOIN ncbi_geneinfo AS i2
  ON i2.entrezID = r2.entrezID AND i2.chr = i1.chr
WHERE
  r1.geneAcc LIKE 'NC%'


/* count pheno/gene associations per pval magnitude */
SELECT
  COUNT(1), pval
FROM (
  SELECT pg.entrezID, FLOOR(MAX(-log10(d.pval))) AS pval
  FROM mcpwas_phenogene AS pg
  JOIN mcpwas_data AS d USING (id)
  WHERE pg.study='MEC'
  GROUP BY pg.entrezID
) AS z
GROUP BY pval
