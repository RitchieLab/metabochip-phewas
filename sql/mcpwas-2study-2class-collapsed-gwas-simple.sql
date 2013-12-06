/* collapsed query (one row per SNP, linked to SNP/pheno GWAS hits) */
SELECT
  qGWAS.SNPID,
  qGWAS.rs,
  qGWAS.chrom,
  qGWAS.pos,
  COUNT(DISTINCT d.id) AS count_associations,
  CONVERT(GROUP_CONCAT(d.study ORDER BY d.pval, d.id) USING 'utf8') AS studies,
  CONVERT(GROUP_CONCAT(d.phenotype ORDER BY d.pval, d.id) USING 'utf8') AS phenotypes,
  CONVERT(GROUP_CONCAT(p.phenotype_long ORDER BY d.pval, d.id) USING 'utf8') AS phenotype_longs,
  CONVERT(GROUP_CONCAT(c.phenotype_class ORDER BY d.pval, d.id) USING 'utf8') AS phenotype_classes,
  CONVERT(GROUP_CONCAT(d.phenotype_transformation ORDER BY d.pval, d.id) USING 'utf8') AS phenotype_transformations,
  CONVERT(GROUP_CONCAT(
    (CASE WHEN ABS(d.pval) >= 1e-2 THEN ROUND(d.pval,2) ELSE CONCAT(
      ROUND(d.pval*POW(10,CEIL(-LOG10(ABS(d.pval)))),2),'e',-CEIL(-LOG10(ABS(d.pval)))
    ) END)
    ORDER BY d.pval, d.id
  ) USING 'utf8') AS pvals,
  CONVERT(GROUP_CONCAT(
    CONCAT(
      (CASE WHEN ABS(d.beta) >= 1e-2 THEN ROUND(d.beta,2) ELSE CONCAT(
        ROUND(d.beta*POW(10,CEIL(-LOG10(ABS(d.beta)))),2),'e',-CEIL(-LOG10(ABS(d.beta)))
      ) END),
      ' (',
      (CASE WHEN ABS(d.SE) >= 1e-2 THEN ROUND(d.SE,2) ELSE CONCAT(
        ROUND(d.SE*POW(10,CEIL(-LOG10(ABS(d.SE)))),2),'e',-CEIL(-LOG10(ABS(d.SE)))
      ) END),
      ')'
    )
    ORDER BY d.pval, d.id
  ) USING 'utf8') AS beta_SEs,
  CONVERT(GROUP_CONCAT(d.N_total ORDER BY d.pval, d.id) USING 'utf8') AS N_totals,
  qGWAS.ref,
  AVG(NULLIF(d.AF_coded_all,0)) AS avg_AFs,
  qGWAS.info,
  qGWAS.ncbi_inside_genes,
  qGWAS.ncbi_upstream_gene,
  qGWAS.ncbi_upstream_gene_dist,
  qGWAS.ncbi_downstream_gene,
  qGWAS.ncbi_downstream_gene_dist,
  qGWAS.gwas_upstream_disease_trait,
  qGWAS.gwas_upstream_distance,
  qGWAS.gwas_downstream_disease_trait,
  qGWAS.gwas_downstream_distance
FROM (
  /* add upstream/downstream GWAS annotation for each SNP */
  SELECT
    qGenes.SNPID,
    qGenes.phenotype_classes,
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
    gcGWASD.disease_trait AS gwas_downstream_disease_trait,
    (gcGWASD.pos - qGenes.pos) AS gwas_downstream_distance
  FROM (
    /* add upstream/inside/downstream gene annotation for each SNP */
    SELECT
      qClasses.SNPID,
      qClasses.phenotype_classes,
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
    FROM (
      /* identify SNPs for which there are multiple classes with significant results in multiple studies each */
      SELECT
        qStudies.SNPID,
        COUNT(DISTINCT qStudies.phenotype_class) AS count_phenotype_classes,
        CONCAT('\t',GROUP_CONCAT(DISTINCT qStudies.phenotype_class ORDER BY qStudies.phenotype_class SEPARATOR '\t'),'\t') AS phenotype_classes
      FROM (
        /* identify SNP-class combinations with significant results in multiple studies */
        SELECT
          dStudies.SNPID,
          cStudies.phenotype_class,
          COUNT(DISTINCT dStudies.study) AS count_studies
        FROM mcpwas_data AS dStudies
        JOIN mcpwas_class AS cStudies USING (phenotype)
        WHERE dStudies.pval > 0
          AND dStudies.pval < 0.01
          AND dStudies.AF_coded_all > 0.01
        GROUP BY dStudies.SNPID, cStudies.phenotype_class
        HAVING count_studies > 1
      ) AS qStudies
      GROUP BY qStudies.SNPID
      HAVING count_phenotype_classes > 0
    ) AS qClasses
    JOIN mcpwas_markers AS mGenes
      ON mGenes.id = qClasses.SNPID
    JOIN mcpwas_rs AS rGenes
      ON rGenes.SNPID = qClasses.SNPID
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
      ON mgGenes.id = qClasses.SNPID
    LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiGenes
      ON ngiGenes.entrezID = mgGenes.entrezID
    GROUP BY qClasses.SNPID
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
  AND INSTR(qGWAS.phenotype_classes, CONCAT('\t',c.phenotype_class,'\t')) > 0
WHERE d.pval > 0 AND d.pval < 0.01
GROUP BY d.SNPID
