/* collapsed query (one row per SNPID) */
SELECT
  m.id,
  CONVERT(CONCAT('rs',r.rs) USING latin1) AS rs,
  m.chrom,
  m.pos,
  COUNT(DISTINCT COALESCE(c.phenotype_superclass,c.phenotype_class)) AS count_phenotype_class,
  GROUP_CONCAT(d.phenotype ORDER BY d.pval, d.id) AS list_phenotype,
  GROUP_CONCAT(p.phenotype_long ORDER BY d.pval, d.id) AS list_phenotype_long,
  GROUP_CONCAT(COALESCE(c.phenotype_superclass,c.phenotype_class) ORDER BY d.pval, d.id) AS list_phenotype_class,
  GROUP_CONCAT(d.study ORDER BY d.pval, d.id) AS list_study,
  GROUP_CONCAT(d.N_total ORDER BY d.pval, d.id) AS list_N_total,
  GROUP_CONCAT(d.pval ORDER BY d.pval, d.id) AS list_pvals,
  GROUP_CONCAT(CONCAT(d.beta,'(',d.SE,')') ORDER BY d.pval, d.id) AS betaSE,
  m.ref AS ref,
  GROUP_CONCAT(d.AF_coded_all ORDER BY d.id) AS AFCodedAll,
  m.info AS info
FROM mcpwas_data AS d
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
  AND (
    c.phenotype_superclass IN ('Lipids', 'ObesityLevel', 'Heme')
    OR c.phenotype_class IN (
      'Activity','Albumin','Alcohol','Arrythmia','ArteryHeartSurgery','AtrialFibrillation','AtrialFlutter','AVBlock',
      'CardiovascularDisease','CHD','Creatinine','CRP','Cystatin','D-Dimer','Diabetes','Diastolic','ECG','Fibrinogen','Glucose','HeartFailure','HeartRate',
      'Height',
      'HormoneUse',
      'Hypertension','Hysterectomy','Insulin','LipidMedications','Menarche','Menopause','MyocardialInfarction','Oophorectomy',
      'Pacemaker','PlateletCount','Pregnancy','PRInterval','QRS','QT','Smoking','Stroke','Systolic','WhiteBloodCount'
    )
  )
WHERE
  d.pval > 0
  AND d.pval < 5e-4
  AND d.AF_coded_all > 0.01
  AND NOT (d.phenotype_transformation = 'None' AND d.effect_size_type = 'Linear regression coefficient')
GROUP BY d.SNPID
HAVING count_phenotype_class > 1
ORDER BY count_phenotype_class DESC
