/* expanded query (one row per data point) for SNPs associated with the same class in multiple studies */
SELECT
  qStudies.marker,
  qStudies.rs,
  qStudies.chrom,
  qStudies.pos,
  qStudies.count_studies,
  ps.study,
  pp.phenotype,
  rc.phenotype_class,
  pt.transform,
  pd.pvalue,
  qStudies.ref,
  rms.AF_coded_all,
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
JOIN mcpwas_permu.data AS pd USING (marker_id)
JOIN mcpwas_permu.phenotype AS pp USING (phenotype_id)
JOIN ritchie_lab.mcpwas_class AS rc USING (phenotype, phenotype_class)
JOIN mcpwas_permu.transform AS pt USING (transform_id)
JOIN mcpwas_permu.study AS ps USING (study_id)
JOIN ritchie_lab.mcpwas_marker_study AS rms
  ON rms.id = qStudies.marker
  AND rms.study = ps.study
  AND rms.substudy = ps.substudy
GROUP BY pd.assoc_id
ORDER BY qStudies.marker, qStudies.phenotype_class, pd.pvalue
