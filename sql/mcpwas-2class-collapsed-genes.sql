/* collapsed query (one row per SNPID) with nearby genes */
SELECT
  m.id,
  CONVERT(CONCAT('rs',r.rs) USING latin1) AS rs,
  m.chrom,
  m.pos,
  qMarker.count_phenotype_class,
  GROUP_CONCAT(d.phenotype ORDER BY d.pval, d.id) AS phenotype,
  GROUP_CONCAT(p.phenotype_long ORDER BY d.pval, d.id) AS phenotype_long,
  GROUP_CONCAT(COALESCE(c.phenotype_superclass,c.phenotype_class) ORDER BY d.pval, d.id) AS phenotype_class,
  GROUP_CONCAT(d.study ORDER BY d.pval, d.id) AS study,
  GROUP_CONCAT(d.phenotype_transformation ORDER BY d.pval, d.id) AS phenotype_transformation,
  GROUP_CONCAT(d.pval ORDER BY d.pval, d.id) AS pval,
  GROUP_CONCAT(CONCAT(d.beta,'(',d.SE,')') ORDER BY d.pval, d.id) AS betaSE,
  GROUP_CONCAT(d.N_total ORDER BY d.pval, d.id) AS N_total,
  m.ref,
  GROUP_CONCAT(d.AF_coded_all ORDER BY d.pval, d.id) AS AFCodedAll,
  m.info,
  qMarker.inside_genes,
  qMarker.upstream_gene,
  qMarker.upstream_gene_dist,
  qMarker.downstream_gene,
  qMarker.downstream_gene_dist
FROM (
  SELECT
    qPheno.SNPID,
    qPheno.count_phenotype_class,
    qPheno.phenotypes,
    qPheno.minPval,
    qPheno.maxPval,
    qPheno.minAFca,
    qPheno.maxAFca,
    GROUP_CONCAT(DISTINCT ngiMarker.symbol ORDER BY ngiMarker.symbol) AS inside_genes,
    ngiMarkerU.symbol AS upstream_gene,
    -(mMarker.pos - MAX(ngrMarkerU.endPos)) AS upstream_gene_dist,
    ngiMarkerD.symbol AS downstream_gene,
    (MIN(ngrMarkerD.startPos) - mMarker.pos) AS downstream_gene_dist
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
      AND dPheno.pval < 5e-4
      AND dPheno.AF_coded_all > 0.01
      AND NOT (dPheno.phenotype_transformation = 'None' AND dPheno.effect_size_type = 'Linear regression coefficient')
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
  AND NOT (d.phenotype_transformation = 'None' AND d.effect_size_type = 'Linear regression coefficient')
JOIN mcpwas_markers AS m
  ON m.id = d.SNPID
JOIN mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
GROUP BY d.SNPID
LIMIT 10
