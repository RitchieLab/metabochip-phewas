/* expanded query (one row per data point, linked to SNP/pheno GWAS hits) */
/* 2:27551325 -> GTF3C2 27548720-27579867 */
SELECT
  m.id,
  (CASE WHEN r.rs IS NULL THEN NULL ELSE CONCAT('rs',r.rs) END) AS rs,
  m.chrom,
  m.pos,
  qMarker.count_phenotype_class,
  d.study,
  d.phenotype,
  p.phenotype_long,
  c.phenotype_class,
  d.phenotype_transformation,
  d.pval,
  d.beta,
  d.SE,
  d.N_total,
  m.ref,
  d.AF_coded_all,
  m.info,
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
    qPheno.SNPID,
    qPheno.count_phenotype_class,
    qPheno.phenotypes,
    qPheno.minPval,
    qPheno.maxPval,
    qPheno.minAFca,
    qPheno.maxAFca,
    GROUP_CONCAT(DISTINCT ngiMarker.symbol ORDER BY ngiMarker.symbol) AS ncbi_inside_genes,
    ngiMarkerU.symbol AS ncbi_upstream_gene,
    -(mMarker.pos - MAX(ngrMarkerU.endPos)) AS ncbi_upstream_gene_dist,
    ngiMarkerD.symbol AS ncbi_downstream_gene,
    (MIN(ngrMarkerD.startPos) - mMarker.pos) AS ncbi_downstream_gene_dist
  FROM (
    SELECT
      dPheno.SNPID,
      COUNT(DISTINCT COALESCE(cPheno.phenotype_superclass,cPheno.phenotype_class)) AS count_phenotype_class,
      CONCAT('\t',GROUP_CONCAT(DISTINCT dPheno.phenotype SEPARATOR '\t'),'\t') AS phenotypes,
      MIN(dPheno.pval) AS minPval,
      MAX(dPheno.pval) AS maxPval,
      MIN(dPheno.AF_coded_all) AS minAFca,
      MAX(dPheno.AF_coded_all) AS maxAFca
    FROM mcpwas_data AS dPheno
    JOIN mcpwas_class AS cPheno
      ON cPheno.phenotype = dPheno.phenotype
      AND (
        cPheno.phenotype_superclass IN ('Lipids', 'ObesityLevel','Heme')
        OR cPheno.phenotype_class IN (
'Activity','Albumin','Alcohol','Arrythmia','ArteryHeartSurgery','AtrialFibrillation','AtrialFlutter','AVBlock',
'CardiovascularDisease','CHD','Creatinine','CRP','Cystatin','D-Dimer','Diabetes','Diastolic','ECG','Fibrinogen','Glucose','HeartFailure','HeartRate',
'Height', 
'HormoneUse','Hypertension','Hysterectomy','Insulin','LipidMedications','Menarche','Menopause','MyocardialInfarction','Oophorectomy',
'Pacemaker','PlateletCount','Pregnancy','PRInterval','QRS','QT','Smoking','Stroke','Systolic','WhiteBloodCount','Wolff-Parkinson-White'
        )
      )
    WHERE dPheno.pval > 0
      AND dPheno.pval < 8e-4
      AND dPheno.AF_coded_all > 0.01
    GROUP BY dPheno.SNPID
    HAVING count_phenotype_class > 1
  ) AS qPheno
  JOIN mcpwas_markers AS mMarker
    ON mMarker.id = qPheno.SNPID
  LEFT JOIN ncbi_geneinfo AS ngiMarkerU
    ON ngiMarkerU.entrezID = mMarker.upstream_entrezID
  LEFT JOIN ncbi_geneinfo AS ngiMarkerD
    ON ngiMarkerD.entrezID = mMarker.downstream_entrezID
  LEFT JOIN ncbi_gene2refseq AS ngrMarkerU
    ON ngrMarkerU.entrezID = ngiMarkerU.entrezID
    AND ngrMarkerU.geneAcc2 = 'NC'
  LEFT JOIN ncbi_gene2refseq AS ngrMarkerD
    ON ngrMarkerD.entrezID = ngiMarkerD.entrezID
    AND ngrMarkerD.geneAcc2 = 'NC'
  LEFT JOIN mcpwas_markergene AS mgMarker
    ON mgMarker.id = qPheno.SNPID
  LEFT JOIN ncbi_geneinfo AS ngiMarker
    ON ngiMarker.entrezID = mgMarker.entrezID
  GROUP BY qPheno.SNPID
) AS qMarker
JOIN mcpwas_data AS d
  ON d.SNPID = qMarker.SNPID
  AND d.pval >= qMarker.minPval
  AND d.pval <= qMarker.maxPval
  AND d.AF_coded_all >= qMarker.minAFca
  AND d.AF_coded_all <= qMarker.maxAFca
  AND INSTR(qMarker.phenotypes, CONCAT('\t',d.phenotype,'\t')) > 0
JOIN mcpwas_markers AS m
  ON m.id = d.SNPID
JOIN mcpwas_rs AS r
  ON r.SNPID = d.SNPID
JOIN mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
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
HAVING gwas_result >= 0
ORDER BY d.SNPID, d.phenotype, d.pval
LIMIT 10

