load data local infile '~/group/projects/MetaboPheWAS/aric_corr3.txt' into table _tmp_corr ignore 1 lines;

load data local infile '~/group/projects/MetaboPheWAS/aric_corr3.txt' into table _tmp_corr fields terminated by '\t' lines terminated by '\n' ignore 1 lines;

select phenotag,substring_index(phenotag,'.',1) as p,phenotype from _tmp_corr having p!=phenotype;

update mcpwas_phenocorr pc join _tmp_corr tc on tc.phenotag=pc.phenotag1 and tc.study=pc.study and tc.substudy=pc.substudy set pc.phenotype1=tc.phenotype;
update mcpwas_phenocorr pc join _tmp_corr tc on tc.phenotag=pc.phenotag2 and tc.study=pc.study and tc.substudy=pc.substudy set pc.phenotype2=tc.phenotype;


SELECT
  sub2.*,
  AVG(d2.N_total) AS avg_N_total_2
FROM (
  SELECT
    sub1.*,
    AVG(d1.N_total) AS avg_N_total_1
  FROM (
    SELECT
      p1.study,
      p1.substudy,
      COALESCE(c1.phenotype_superclass,c1.phenotype_class) AS phenotype_class_1,
      pc.phenotype1 AS phenotype_1,
      p1.phenotype_long AS phenotype_long_1,
      pc.transform1 AS phenotype_transformation_1,
      COALESCE(c2.phenotype_superclass,c2.phenotype_class) AS phenotype_class_2,
      pc.phenotype2 AS phenotype_2,
      p2.phenotype_long AS phenotype_long_2,
      pc.transform2 AS phenotype_transformation_2,
      pc.correlation
    FROM mcpwas_pheno AS p1
    JOIN mcpwas_class AS c1
      ON c1.phenotype = p1.phenotype
    JOIN mcpwas_phenocorr AS pc
      ON pc.study = p1.study AND pc.substudy = p1.substudy AND pc.phenotype1 = c1.phenotype
    LEFT JOIN mcpwas_pheno AS p2
      ON p2.study = p1.study AND p2.substudy = p1.substudy AND p2.phenotype = pc.phenotype2
    LEFT JOIN mcpwas_class AS c2
      ON c2.phenotype = p2.phenotype
    WHERE ABS(pc.correlation) >= 0.6
    HAVING phenotype_class_1 != phenotype_class_2
  ) AS sub1
  JOIN mcpwas_data AS d1
    ON d1.study = sub1.study
    AND d1.substudy = sub1.substudy
    AND d1.phenotype = sub1.phenotype_1
  GROUP BY sub1.study, sub1.substudy, sub1.phenotype_1, sub1.phenotype_2
) AS sub2
JOIN mcpwas_data AS d2
  ON d2.study = sub2.study
  AND d2.substudy = sub2.substudy
  AND d2.phenotype = sub2.phenotype_2
GROUP BY sub2.study, sub2.substudy, sub2.phenotype_1, sub2.phenotype_2


SELECT
  p1.study,
  p1.substudy,
  COALESCE(c1.phenotype_superclass,c1.phenotype_class) AS phenotype_class_1,
  pc.phenotype1 AS phenotype_1,
  p1.phenotype_long AS phenotype_long_1,
  pc.transform1 AS phenotype_transformation_1,
  (
    SELECT AVG(N_total) FROM mcpwas_data AS d1
    ON d1.study = p1.study AND d1.substudy = p1.substudy AND d1.phenotype = p1.phenotype
  ) AS avg_N_total_1,
  COALESCE(c2.phenotype_superclass,c2.phenotype_class) AS phenotype_class_2,
  pc.phenotype2 AS phenotype_2,
  p2.phenotype_long AS phenotype_long_2,
  pc.transform2 AS phenotype_transformation_2,
  (
    SELECT AVG(N_total) FROM mcpwas_data AS d2
    ON d2.study = p2.study AND d2.substudy = p2.substudy AND d2.phenotype = p2.phenotype
  ) AS avg_N_total_2,
  pc.correlation
FROM mcpwas_pheno AS p1
JOIN mcpwas_class AS c1
  ON c1.phenotype = p1.phenotype
JOIN mcpwas_phenocorr AS pc
  ON pc.study = p1.study AND pc.substudy = p1.substudy AND pc.phenotype1 = c1.phenotype
LEFT JOIN mcpwas_pheno AS p2
  ON p2.study = p1.study AND p2.substudy = p1.substudy AND p2.phenotype = pc.phenotype2
LEFT JOIN mcpwas_class AS c2
  ON c2.phenotype = p2.phenotype
WHERE ABS(pc.correlation) >= 0.6
HAVING phenotype_class_1 != phenotype_class_2


SELECT
  p1.study,
  p1.substudy,
  COALESCE(c1.phenotype_superclass,c1.phenotype_class) AS phenotype_class_1,
  pc.phenotype1 AS phenotype_1,
  p1.phenotype_long AS phenotype_long_1,
  pc.transform1 AS phenotype_transformation_1,
  p1d.avg_N_total AS avg_N_total_1,
  COALESCE(c2.phenotype_superclass,c2.phenotype_class) AS phenotype_class_2,
  pc.phenotype2 AS phenotype_2,
  p2.phenotype_long AS phenotype_long_2,
  pc.transform2 AS phenotype_transformation_2,
  p2d.avg_N_total AS avg_N_total_2,
  pc.correlation
FROM mcpwas_pheno AS p1
JOIN mcpwas_class AS c1
  ON c1.phenotype = p1.phenotype
JOIN (
  SELECT study,substudy,phenotype,AVG(N_total) AS avg_N_total
  FROM mcpwas_data GROUP BY study,substudy,phenotype
) AS p1d
  ON p1d.study = p1.study AND p1d.substudy = p1.substudy AND p1d.phenotype = p1.phenotype
JOIN mcpwas_phenocorr AS pc
  ON pc.study = p1.study AND pc.substudy = p1.substudy AND pc.phenotype1 = c1.phenotype
  AND ABS(pc.correlation) >= 0.6
LEFT JOIN mcpwas_pheno AS p2
  ON p2.study = p1.study AND p2.substudy = p1.substudy AND p2.phenotype = pc.phenotype2
LEFT JOIN mcpwas_class AS c2
  ON c2.phenotype = p2.phenotype
JOIN (
  SELECT study,substudy,phenotype,AVG(N_total) AS avg_N_total
  FROM mcpwas_data GROUP BY study,substudy,phenotype
) AS p2d
  ON p2d.study = p2.study AND p2d.substudy = p2.substudy AND p2d.phenotype = p2.phenotype
HAVING phenotype_class_1 != phenotype_class_2

