/* collapsed query (one row per SNP point, linked to SNP/pheno GWAS hits) */
SELECT
  qGenes.SNPID,
  qGenes.rs,
  qGenes.chrom,
  qGenes.pos,
  qGenes.count_phenotype_class,
  GROUP_CONCAT(
    d.study
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS study,
  GROUP_CONCAT(
    d.phenotype
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS phenotype,
  GROUP_CONCAT(
    p.phenotype_long
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS phenotype_long,
  GROUP_CONCAT(
    c.phenotype_class
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS phenotype_class,
  GROUP_CONCAT(
    d.phenotype_transformation
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS phenotype_transformation,
  GROUP_CONCAT(
    d.pval
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS pval,
  GROUP_CONCAT(
    d.beta
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS beta,
  GROUP_CONCAT(
    d.SE
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS SE,
  GROUP_CONCAT(
    d.N_total
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS N_total,
  qGenes.ref,
  GROUP_CONCAT(
    d.AF_coded_all
    ORDER BY LEAST(g.upstream_gwas_hit_distance,g.downstream_gwas_hit_distance), d.id, g.id
  ) AS AF_coded_all,
  qGenes.info,
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
  qGenes.ncbi_inside_genes,
  qGenes.ncbi_upstream_gene,
  qGenes.ncbi_upstream_gene_dist,
  qGenes.ncbi_downstream_gene,
  qGenes.ncbi_downstream_gene_dist
FROM (
    /* add upstream/inside/downstream gene annotation for each SNP */
    SELECT
      qClasses.SNPID,
      qClasses.count_phenotype_class,
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
        COUNT(DISTINCT qStudies.phenotype_class) AS count_phenotype_class,
        CONCAT('\t',GROUP_CONCAT(DISTINCT qStudies.phenotype_class ORDER BY qStudies.phenotype_class SEPARATOR '\t'),'\t') AS phenotype_classes,
        CONCAT('\t',GROUP_CONCAT(DISTINCT pgClasses.disease_trait SEPARATOR '\t'),'\t') AS class_disease_traits
      FROM (
        /* identify SNP-class combinations with significant results in multiple studies */
        SELECT
          dStudies.SNPID,
          cStudies.phenotype_class,
          COUNT(DISTINCT dStudies.study) AS count_study
        FROM mcpwas_data AS dStudies
        JOIN mcpwas_class AS cStudies USING (phenotype)
        WHERE dStudies.pval > 0
          AND dStudies.pval < 0.01
          AND dStudies.AF_coded_all > 0.01
        GROUP BY dStudies.SNPID, cStudies.phenotype_class
        HAVING count_study > 1
      ) AS qStudies
      LEFT JOIN mcpwas_phenogwas AS pgClasses
        ON pgClasses.phenotype = qStudies.phenotype_class
      GROUP BY qStudies.SNPID
      HAVING count_phenotype_class > 1
    ) AS qClasses
    JOIN mcpwas_markers AS mGenes
      ON mGenes.id = qClasses.SNPID
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
    GROUP BY qClasses.SNPID
) AS qGenes
JOIN mcpwas_data AS d
  ON d.SNPID = qGenes.SNPID
JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
  AND INSTR(qGenes.phenotype_classes, CONCAT('\t',c.phenotype_class,'\t')) > 0
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
GROUP BY d.SNPID
ORDER BY d.SNPID, c.phenotype_class, d.pval
