/* expanded query (one row per data point, linked to SNP/pheno GWAS hits) */
SELECT
  qGWAS.SNPID,
  qGWAS.rs,
  qGWAS.chrom,
  qGWAS.pos,
  d.study,
  d.phenotype,
  p.phenotype_long,
  c.phenotype_class,
  d.phenotype_transformation,
  CONVERT((CASE WHEN ABS(d.pval) >= 1e-2 THEN ROUND(d.pval,2) ELSE CONCAT(
    ROUND(d.pval*POW(10,CEIL(-LOG10(ABS(d.pval)))),2),'e',-CEIL(-LOG10(ABS(d.pval)))
  ) END) USING 'utf8') AS pval,
  CONVERT((CASE WHEN ABS(d.beta) >= 1e-2 THEN ROUND(d.beta,2) ELSE CONCAT(
    ROUND(d.beta*POW(10,CEIL(-LOG10(ABS(d.beta)))),2),'e',-CEIL(-LOG10(ABS(d.beta)))
  ) END) USING 'utf8') AS beta,
  CONVERT((CASE WHEN ABS(d.SE) >= 1e-2 THEN ROUND(d.SE,2) ELSE CONCAT(
    ROUND(d.SE*POW(10,CEIL(-LOG10(ABS(d.SE)))),2),'e',-CEIL(-LOG10(ABS(d.SE)))
  ) END) USING 'utf8') AS SE,
  d.N_total,
  qGWAS.ref,
  d.AF_coded_all,
  qGWAS.info,
  qGWAS.ncbi_inside_genes,
  qGWAS.ncbi_upstream_gene,
  qGWAS.ncbi_upstream_gene_dist,
  qGWAS.ncbi_downstream_gene,
  qGWAS.ncbi_downstream_gene_dist,
  qGWAS.gwas_upstream_disease_trait,
  qGWAS.gwas_upstream_distance,
  qGWAS.gwas_upstream_pubmedid,
  qGWAS.gwas_downstream_disease_trait,
  qGWAS.gwas_downstream_distance,
  qGWAS.gwas_downstream_pubmedid
FROM (
  /* add upstream/downstream GWAS annotation for each SNP */
  SELECT
    qGenes.SNPID,
    qGenes.rs,
    qGenes.chrom,
    qGenes.pos,
    qGenes.ref,
    qGenes.info,
    qGenes.ncbi_inside_genes,
    qGenes.ncbi_upstream_gene,
    qGenes.ncbi_upstream_gene_dist,
    qGenes.ncbi_downstream_gene,
    qGenes.ncbi_downstream_gene_dist,
    gcGWASU.disease_trait AS gwas_upstream_disease_trait,
    -(gcGWASU.pos - qGenes.pos) AS gwas_upstream_distance,
    gcGWASU.pubmedid AS gwas_upstream_pubmedid,
    gcGWASD.disease_trait AS gwas_downstream_disease_trait,
    (gcGWASD.pos - qGenes.pos) AS gwas_downstream_distance,
    gcGWASD.pubmedid AS gwas_downstream_pubmedid
  FROM (
    /* add upstream/inside/downstream gene annotation for each SNP */
    SELECT
      mGenes.id AS SNPID,
      CONVERT(CASE WHEN rGenes.rs IS NULL THEN NULL ELSE CONCAT('rs',rGenes.rs) END USING utf8) AS rs,
      mGenes.chrom,
      mGenes.pos,
      mGenes.ref,
      mGenes.info,
      GROUP_CONCAT(DISTINCT ngiGenes.symbol ORDER BY ngiGenes.symbol) AS ncbi_inside_genes,
      ngiGenesU.symbol AS ncbi_upstream_gene,
      -(mGenes.pos - MAX(ngrGenesU.endPos)) AS ncbi_upstream_gene_dist,
      ngiGenesD.symbol AS ncbi_downstream_gene,
      (MIN(ngrGenesD.startPos) - mGenes.pos) AS ncbi_downstream_gene_dist,
      (SELECT MAX(gcGenesU.pos) FROM gwas_catalog AS gcGenesU WHERE gcGenesU.chr = mGenes.chrom AND gcGenesU.pos <= mGenes.pos AND gcGenesU.pos >= mGenes.pos - 100000) AS gwas_upstream_pos,
      (SELECT MIN(gcGenesD.pos) FROM gwas_catalog AS gcGenesD WHERE gcGenesD.chr = mGenes.chrom AND gcGenesD.pos >= mGenes.pos AND gcGenesD.pos <= mGenes.pos + 100000) AS gwas_downstream_pos
    FROM mcpwas_markers AS mGenes
    JOIN mcpwas_rs AS rGenes
      ON rGenes.SNPID = mGenes.id
    LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiGenesU
      ON ngiGenesU.entrezID = mGenes.upstream_entrezID
    LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiGenesD
      ON ngiGenesD.entrezID = mGenes.downstream_entrezID
    LEFT JOIN ritchie_lab.ncbi_gene2refseq AS ngrGenesU
      ON ngrGenesU.entrezID = ngiGenesU.entrezID
      AND ngrGenesU.geneAcc2 = 'NC'
    LEFT JOIN ritchie_lab.ncbi_gene2refseq AS ngrGenesD
      ON ngrGenesD.entrezID = ngiGenesD.entrezID
      AND ngrGenesD.geneAcc2 = 'NC'
    LEFT JOIN mcpwas_markergene AS mgGenes
      ON mgGenes.id = mGenes.id
    LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiGenes
      ON ngiGenes.entrezID = mgGenes.entrezID
    GROUP BY mGenes.id
  ) AS qGenes
  LEFT JOIN gwas_catalog AS gcGWASU
    ON gcGWASU.chr = qGenes.chrom
    AND gcGWASU.pos = qGenes.gwas_upstream_pos
  LEFT JOIN gwas_catalog AS gcGWASD
    ON gcGWASD.chr = qGenes.chrom
    AND gcGWASD.pos = qGenes.gwas_downstream_pos
  GROUP BY qGenes.SNPID
) AS qGWAS
JOIN mcpwas_data AS d
  ON d.SNPID = qGWAS.SNPID
JOIN mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
WHERE d.pval > 0
  AND d.pval < 1e-4
  AND d.AF_coded_all > 0.01
GROUP BY d.id
ORDER BY d.SNPID, c.phenotype_class, d.pval
