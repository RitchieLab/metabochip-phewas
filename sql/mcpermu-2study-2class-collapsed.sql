/* expanded query (one row per data point) for SNPs associated with the same class in multiple studies for multiple classes */
SELECT
  qClasses.marker,
  qClasses.rs,
  qClasses.chrom,
  qClasses.pos,
  qClasses.count_phenotype_classes,
  COUNT(1) AS count_associations,
  CONVERT(GROUP_CONCAT(ps.study ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  CONVERT(GROUP_CONCAT(pp.phenotype ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  CONVERT(GROUP_CONCAT(rc.phenotype_class ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  CONVERT(GROUP_CONCAT(pt.transform ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  CONVERT(GROUP_CONCAT(pd.pvalue ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  qClasses.ref,
  CONVERT(GROUP_CONCAT(rms.AF_coded_all ORDER BY rc.phenotype_class, pd.pvalue, pd.assoc_id) USING 'utf8'),
  qClasses.info
FROM (
  SELECT
    qStudies.marker_id,
    qStudies.marker,
    COUNT(DISTINCT qStudies.phenotype_class) AS count_phenotype_classes,
    CONCAT('\t',GROUP_CONCAT(DISTINCT qStudies.phenotype_class SEPARATOR '\t'),'\t') AS phenotype_classes,
    qStudies.rs,
    qStudies.chrom,
    qStudies.pos,
    qStudies.ref,
    qStudies.info
  FROM (
    SELECT
      pdStudies.marker_id,
      pmStudies.marker,
      CONVERT(CASE WHEN rrStudies.rs IS NULL THEN NULL ELSE CONCAT('rs',rrStudies.rs) END USING utf8) AS rs,
      rmStudies.chrom,
      rmStudies.pos,
      COUNT(DISTINCT pdStudies.study_id) AS count_studies,
      rcStudies.phenotype_class,
      rmStudies.ref,
      rmStudies.info
    FROM mcpwas_permu.data AS pdStudies
    JOIN mcpwas_permu.marker AS pmStudies USING (marker_id)
    JOIN mcpwas_permu.study AS psStudies USING (study_id)
    JOIN mcpwas_permu.phenotype AS ppStudies USING (phenotype_id)
    JOIN ritchie_lab.mcpwas_markers AS rmStudies
      ON rmStudies.id = pmStudies.marker
    JOIN ritchie_lab.mcpwas_rs AS rrStudies
      ON rrStudies.SNPID = pmStudies.marker
    JOIN ritchie_lab.mcpwas_marker_study AS rmsStudies
      ON rmsStudies.id = pmStudies.marker
      AND rmsStudies.study = psStudies.study
      AND rmsStudies.substudy = psStudies.substudy
    JOIN ritchie_lab.mcpwas_class AS rcStudies USING (phenotype)
    WHERE pdStudies.permu_id = 1
      AND pdStudies.pvalue > 0
      AND pdStudies.pvalue < 1e-4
      AND rmsStudies.AF_coded_all > 0.01
      AND rcStudies.phenotype_class IS NOT NULL
    GROUP BY pdStudies.marker_id, rcStudies.phenotype_class
    HAVING count_studies > 1
  ) AS qStudies
  GROUP BY qStudies.marker_id
  HAVING count_phenotype_classes > 1
) AS qClasses
JOIN mcpwas_permu.data AS pd USING (marker_id)
JOIN mcpwas_permu.phenotype AS pp USING (phenotype_id)
JOIN ritchie_lab.mcpwas_class AS rc
  ON rc.phenotype = pp.phenotype
  AND rc.phenotype_class IS NOT NULL
  AND INSTR(qClasses.phenotype_classes, CONCAT('\t',rc.phenotype_class,'\t')) > 0
JOIN mcpwas_permu.transform AS pt USING (transform_id)
JOIN mcpwas_permu.study AS ps USING (study_id)
JOIN ritchie_lab.mcpwas_marker_study AS rms
  ON rms.id = qClasses.marker
  AND rms.study = ps.study
  AND rms.substudy = ps.substudy
GROUP BY qClasses.marker_id
ORDER BY qClasses.marker
