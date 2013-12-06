/* expanded query (one row per data point, linked to SNP/pheno GWAS hits) */
SELECT
  qMarker.SNPID,
  qMarker.rs,
  qMarker.chrom,
  qMarker.pos,
  qMarker.count_studies,
  d.study,
  d.phenotype,
  p.phenotype_long,
  c.phenotype_class,
  d.phenotype_transformation,
  d.pval,
  d.beta,
  d.SE,
  d.N_total,
  qMarker.ref,
  d.AF_coded_all,
  qMarker.info,
  MAX(CASE
    WHEN pg.disease_trait IS NULL THEN 0
    WHEN g.id IS NULL THEN 1
    WHEN g.upstream_gwas_hit_distance = 0 AND g.downstream_gwas_hit_distance = 0 THEN 3
    ELSE 2
  END) AS gwas_result,
  GROUP_CONCAT(
    (CASE
      WHEN pg.disease_trait IS NULL THEN NULL
      WHEN g.id IS NULL THEN NULL
      ELSE pg.disease_trait
    END)
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS gwas_disease_traits,
  GROUP_CONCAT(
    (CASE
      WHEN pg.disease_trait IS NULL THEN NULL
      WHEN g.id IS NULL THEN NULL
      WHEN (
        g.upstream_gwas_hit_distance < 10000
        AND pg.disease_trait IN (
          g.upstream_gwas_hit_disease_trait_1,
          g.upstream_gwas_hit_disease_trait_2,
          g.upstream_gwas_hit_disease_trait_3,
          g.upstream_gwas_hit_disease_trait_4,
          g.upstream_gwas_hit_disease_trait_5,
          g.upstream_gwas_hit_disease_trait_6,
          g.upstream_gwas_hit_disease_trait_7
        )
      ) THEN -g.upstream_gwas_hit_distance
      ELSE g.downstream_gwas_hit_distance
    END)
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS gwas_hit_distances,
  GROUP_CONCAT(
    g.upstream_closest_gene
    ORDER BY g.upstream_closest_gene_distance
  ) AS gwas_upstream_genes,
  GROUP_CONCAT(
    -g.upstream_closest_gene_distance
    ORDER BY g.upstream_closest_gene_distance
  ) AS gwas_upstream_gene_distances,
  GROUP_CONCAT(
    g.downstream_closest_gene
    ORDER BY g.downstream_closest_gene_distance
  ) AS gwas_downstream_genes,
  GROUP_CONCAT(
    g.downstream_closest_gene_distance
    ORDER BY g.downstream_closest_gene_distance
  ) AS gwas_downstream_gene_distances,
  qMarker.ncbi_inside_genes,
  qMarker.ncbi_upstream_gene,
  qMarker.ncbi_upstream_gene_dist,
  qMarker.ncbi_downstream_gene,
  qMarker.ncbi_downstream_gene_dist
FROM (
  SELECT
    qSNP.SNPID,
    qSNP.phenotype_class,
    qSNP.count_studies,
    CONVERT(CASE WHEN rMarker.rs IS NULL THEN NULL ELSE CONCAT('rs',rMarker.rs) END USING utf8) AS rs,
    mMarker.chrom,
    mMarker.pos,
    mMarker.ref,
    mMarker.info,
    GROUP_CONCAT(DISTINCT ngiMarker.symbol ORDER BY ngiMarker.symbol) AS ncbi_inside_genes,
    ngiMarkerU.symbol AS ncbi_upstream_gene,
    -(mMarker.pos - MAX(ngrMarkerU.endPos)) AS ncbi_upstream_gene_dist,
    ngiMarkerD.symbol AS ncbi_downstream_gene,
    (MIN(ngrMarkerD.startPos) - mMarker.pos) AS ncbi_downstream_gene_dist
  FROM (
    SELECT
      dSNP.SNPID,
      cSNP.phenotype_class,
      COUNT(DISTINCT dSNP.study) AS count_studies
    FROM mcpwas_data AS dSNP
    JOIN mcpwas_class AS cSNP USING (phenotype)
    WHERE dSNP.pval > 0
      AND dSNP.pval < 1e-4
      AND dSNP.AF_coded_all > 0.01
    GROUP BY dSNP.SNPID, cSNP.phenotype_class
    HAVING count_studies > 1
  ) AS qSNP
  JOIN mcpwas_markers AS mMarker
    ON mMarker.id = qSNP.SNPID
  JOIN mcpwas_rs AS rMarker
    ON rMarker.SNPID = qSNP.SNPID
  LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiMarkerU
    ON ngiMarkerU.entrezID = mMarker.upstream_entrezID
  LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiMarkerD
    ON ngiMarkerD.entrezID = mMarker.downstream_entrezID
  LEFT JOIN ritchie_lab.ncbi_gene2refseq AS ngrMarkerU
    ON ngrMarkerU.entrezID = ngiMarkerU.entrezID
    AND ngrMarkerU.geneAcc2 = 'NC'
  LEFT JOIN ritchie_lab.ncbi_gene2refseq AS ngrMarkerD
    ON ngrMarkerD.entrezID = ngiMarkerD.entrezID
    AND ngrMarkerD.geneAcc2 = 'NC'
  LEFT JOIN mcpwas_markergene AS mgMarker
    ON mgMarker.id = qSNP.SNPID
  LEFT JOIN ritchie_lab.ncbi_geneinfo AS ngiMarker
    ON ngiMarker.entrezID = mgMarker.entrezID
  GROUP BY qSNP.SNPID, qSNP.phenotype_class
) AS qMarker
JOIN mcpwas_data AS d
  ON d.SNPID = qMarker.SNPID
JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
  AND c.phenotype_class = qMarker.phenotype_class
JOIN mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
LEFT JOIN mcpwas_phenogwas AS pg
  ON pg.phenotype = c.phenotype_class
LEFT JOIN mcpwas_gwas AS g
  ON g.id = d.SNPID
  AND pg.disease_trait IS NOT NULL
  AND (
    (
      g.upstream_gwas_hit_distance < 10000
      AND pg.disease_trait IN (
        g.upstream_gwas_hit_disease_trait_1,
        g.upstream_gwas_hit_disease_trait_2,
        g.upstream_gwas_hit_disease_trait_3,
        g.upstream_gwas_hit_disease_trait_4,
        g.upstream_gwas_hit_disease_trait_5,
        g.upstream_gwas_hit_disease_trait_6,
        g.upstream_gwas_hit_disease_trait_7
      )
    ) OR (
      g.downstream_gwas_hit_distance < 10000
      AND pg.disease_trait IN (
        g.downstream_gwas_hit_disease_trait_1,
        g.downstream_gwas_hit_disease_trait_2,
        g.downstream_gwas_hit_disease_trait_3,
        g.downstream_gwas_hit_disease_trait_4,
        g.downstream_gwas_hit_disease_trait_5,
        g.downstream_gwas_hit_disease_trait_6,
        g.downstream_gwas_hit_disease_trait_7
      )
    )
  )
GROUP BY d.id
ORDER BY d.SNPID, c.phenotype_class, d.pval
