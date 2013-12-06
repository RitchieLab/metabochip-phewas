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
  GROUP_CONCAT(
    (CASE WHEN d.pval >= 1e-2 THEN ROUND(d.pval,2) ELSE CONCAT(
      ROUND(d.pval*POW(10,CEIL(-LOG10(d.pval))),2),'e',-CEIL(-LOG10(d.pval))
    ) END)
    ORDER BY d.pval, d.id
  ) AS pvals,
  GROUP_CONCAT(
    CONCAT(
      (CASE WHEN d.beta >= 1e-2 THEN ROUND(d.beta,2) ELSE CONCAT(
        ROUND(d.beta*POW(10,CEIL(-LOG10(d.beta))),2),'e',-CEIL(-LOG10(d.beta))
      ) END),
      ' (',
      (CASE WHEN d.SE >= 1e-2 THEN ROUND(d.SE,2) ELSE CONCAT(
        ROUND(d.SE*POW(10,CEIL(-LOG10(d.SE))),2),'e',-CEIL(-LOG10(d.SE))
      ) END),
      ')'
    )
    ORDER BY d.pval, d.id
  ) AS beta_SEs,
  CONVERT(GROUP_CONCAT(d.N_total ORDER BY d.pval, d.id) USING 'utf8') AS N_totals,
  qGWAS.ref,
  AVG(NULLIF(d.AF_coded_all,0)) AS avg_AFs,
  qGWAS.info,
  qGWAS.ncbi_inside_genes,
  qGWAS.ncbi_upstream_gene,
  qGWAS.ncbi_upstream_gene_dist,
  qGWAS.ncbi_downstream_gene,
  qGWAS.ncbi_downstream_gene_dist,
  qGWAS.gwas_result,
  REPLACE(qGWAS.class_disease_traits,'\t','/') AS class_disease_traits,
  qGWAS.gwas_upstream_disease_traits,
  qGWAS.gwas_upstream_hit_distances,
  qGWAS.gwas_downstream_disease_traits,
  qGWAS.gwas_downstream_hit_distances,
  qGWAS.gwas_upstream_genes,
  qGWAS.gwas_upstream_gene_distances,
  qGWAS.gwas_downstream_genes,
  qGWAS.gwas_downstream_gene_distances
FROM (
  /* add previous GWAS hit annotation for each SNP/class */
  SELECT
    qGenes.SNPID,
    qGenes.phenotype_classes,
    qGenes.class_disease_traits,
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
    (CASE MAX(CASE
        WHEN qGenes.class_disease_traits IS NULL THEN 0
        WHEN gGWAS.id IS NULL THEN 1
        WHEN gGWAS.upstream_gwas_hit_distance = 0 AND gGWAS.downstream_gwas_hit_distance = 0 THEN 3
        ELSE 2
      END)
      WHEN 0 THEN 'no disease trait'
      WHEN 1 THEN 'no nearby GWAS hits'
      WHEN 2 THEN 'nearby GWAS hit(s)'
      WHEN 3 THEN 'exact GWAS hit(s)'
      ELSE 'error!'
    END) AS gwas_result,
    GROUP_CONCAT(DISTINCT
      CONCAT_WS('/',
        gGWAS.upstream_gwas_hit_disease_trait_1,
        gGWAS.upstream_gwas_hit_disease_trait_2,
        gGWAS.upstream_gwas_hit_disease_trait_3,
        gGWAS.upstream_gwas_hit_disease_trait_4,
        gGWAS.upstream_gwas_hit_disease_trait_5,
        gGWAS.upstream_gwas_hit_disease_trait_6,
        gGWAS.upstream_gwas_hit_disease_trait_7
      )
      ORDER BY -gGWAS.upstream_gwas_hit_distance
    ) AS gwas_upstream_disease_traits,
    GROUP_CONCAT(
      DISTINCT -gGWAS.upstream_gwas_hit_distance
      ORDER BY -gGWAS.upstream_gwas_hit_distance
    ) AS gwas_upstream_hit_distances,
    GROUP_CONCAT(DISTINCT
      CONCAT_WS('/',
        gGWAS.downstream_gwas_hit_disease_trait_1,
        gGWAS.downstream_gwas_hit_disease_trait_2,
        gGWAS.downstream_gwas_hit_disease_trait_3,
        gGWAS.downstream_gwas_hit_disease_trait_4,
        gGWAS.downstream_gwas_hit_disease_trait_5,
        gGWAS.downstream_gwas_hit_disease_trait_6,
        gGWAS.downstream_gwas_hit_disease_trait_7
      )
      ORDER BY gGWAS.downstream_gwas_hit_distance
    ) AS gwas_downstream_disease_traits,
    GROUP_CONCAT(
      DISTINCT gGWAS.downstream_gwas_hit_distance
      ORDER BY gGWAS.downstream_gwas_hit_distance
    ) AS gwas_downstream_hit_distances,
    GROUP_CONCAT(DISTINCT gGWAS.upstream_closest_gene ORDER BY gGWAS.upstream_closest_gene_distance) AS gwas_upstream_genes,
    GROUP_CONCAT(DISTINCT -gGWAS.upstream_closest_gene_distance ORDER BY gGWAS.upstream_closest_gene_distance) AS gwas_upstream_gene_distances,
    GROUP_CONCAT(DISTINCT gGWAS.downstream_closest_gene ORDER BY gGWAS.downstream_closest_gene_distance) AS gwas_downstream_genes,
    GROUP_CONCAT(DISTINCT gGWAS.downstream_closest_gene_distance ORDER BY gGWAS.downstream_closest_gene_distance) AS gwas_downstream_gene_distances
  FROM (
    /* add upstream/inside/downstream gene annotation for each SNP */
    SELECT
      qClasses.SNPID,
      qClasses.phenotype_classes,
      qClasses.class_disease_traits,
      CONVERT(CASE WHEN rGenes.rs IS NULL THEN NULL ELSE CONCAT('rs',rGenes.rs) END USING utf8) AS rs,
      mGenes.chrom,
      mGenes.pos,
      mGenes.ref,
      mGenes.info,
      GROUP_CONCAT(DISTINCT ngiGenes.symbol ORDER BY ngiGenes.symbol) AS ncbi_inside_genes,
      ngiGenesU.symbol AS ncbi_upstream_gene,
      -(mGenes.pos - MAX(ngrGenesU.endPos)) AS ncbi_upstream_gene_dist,
      ngiGenesD.symbol AS ncbi_downstream_gene,
      (MIN(ngrGenesD.startPos) - mGenes.pos) AS ncbi_downstream_gene_dist
    FROM (
      /* identify SNPs for which there are multiple classes with significant results in multiple studies each */
      SELECT
        qStudies.SNPID,
        COUNT(DISTINCT qStudies.phenotype_class) AS count_phenotype_classes,
        CONCAT('\t',GROUP_CONCAT(DISTINCT qStudies.phenotype_class ORDER BY qStudies.phenotype_class SEPARATOR '\t'),'\t') AS phenotype_classes,
        CONCAT('\t',GROUP_CONCAT(DISTINCT pgClasses.disease_trait SEPARATOR '\t'),'\t') AS class_disease_traits
      FROM (
        /* identify SNP-class combinations with significant results in multiple studies */
        SELECT
          dStudies.SNPID,
          cStudies.phenotype_class,
          COUNT(DISTINCT dStudies.study) AS count_studies
        FROM ritchie_metabopwas.mcpwas_data AS dStudies
        JOIN ritchie_metabopwas.mcpwas_class AS cStudies USING (phenotype)
        WHERE dStudies.pval > 0
          AND dStudies.pval < 0.01
          AND dStudies.AF_coded_all > 0.01
        GROUP BY dStudies.SNPID, cStudies.phenotype_class
        HAVING count_studies > 1
      ) AS qStudies
      LEFT JOIN ritchie_metabopwas.mcpwas_phenogwas AS pgClasses
        ON pgClasses.phenotype = qStudies.phenotype_class
      GROUP BY qStudies.SNPID
      HAVING count_phenotype_classes >= 0
    ) AS qClasses
    JOIN ritchie_metabopwas.mcpwas_markers AS mGenes
      ON mGenes.id = qClasses.SNPID
    JOIN ritchie_metabopwas.mcpwas_rs AS rGenes
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
    LEFT JOIN ritchie_metabopwas.mcpwas_markergene AS mgGenes
      ON mgGenes.id = qClasses.SNPID
    LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiGenes
      ON ngiGenes.entrezID = mgGenes.entrezID
    GROUP BY qClasses.SNPID
  ) AS qGenes
  LEFT JOIN ritchie_metabopwas.mcpwas_gwas AS gGWAS
    ON gGWAS.id = qGenes.SNPID
    AND (
      (
        gGWAS.upstream_gwas_hit_distance < 10000
        AND (
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_1,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_2,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_3,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_4,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_5,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_6,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.upstream_gwas_hit_disease_trait_7,'\t')) > 0)
        )
      ) OR (
        gGWAS.downstream_gwas_hit_distance < 10000
        AND (
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_1,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_2,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_3,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_4,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_5,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_6,'\t')) > 0) OR
          (INSTR(qGenes.class_disease_traits, CONCAT('\t',gGWAS.downstream_gwas_hit_disease_trait_7,'\t')) > 0)
        )
      )
    )
  GROUP BY qGenes.SNPID
) AS qGWAS
JOIN ritchie_metabopwas.mcpwas_data AS d
  ON d.SNPID = qGWAS.SNPID
JOIN ritchie_metabopwas.mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
JOIN ritchie_metabopwas.mcpwas_class AS c
  ON c.phenotype = d.phenotype
  AND INSTR(qGWAS.phenotype_classes, CONCAT('\t',c.phenotype_class,'\t')) > 0
GROUP BY d.SNPID
