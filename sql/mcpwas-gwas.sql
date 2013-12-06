SELECT
  m.id,
  m.chrom,
  m.pos,
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
  g.upstream_gwas_hit,
  g.upstream_gwas_hit_pvalue,
  g.upstream_gwas_hit_distance,
  g.upstream_gwas_hit_disease_trait_1,
  g.upstream_gwas_hit_disease_trait_2,
  g.upstream_gwas_hit_disease_trait_3,
  g.upstream_gwas_hit_disease_trait_4,
  g.upstream_gwas_hit_disease_trait_5,
  g.upstream_gwas_hit_disease_trait_6,
  g.upstream_gwas_hit_disease_trait_7,
  g.downstream_gwas_hit,
  g.downstream_gwas_hit_pvalue,
  g.downstream_gwas_hit_distance,
  g.downstream_gwas_hit_disease_trait_1,
  g.downstream_gwas_hit_disease_trait_2,
  g.downstream_gwas_hit_disease_trait_3,
  g.downstream_gwas_hit_disease_trait_4,
  g.downstream_gwas_hit_disease_trait_5,
  g.downstream_gwas_hit_disease_trait_6,
  g.downstream_gwas_hit_disease_trait_7,
  g.gene,
  g.upstream_closest_gene,
  g.upstream_closest_gene_distance,
  g.downstream_closest_gene,
  g.downstream_closest_gene_distance
FROM mcpwas_data AS d
LEFT JOIN mcpwas_markers AS m
  ON m.id = d.SNPID
LEFT JOIN mcpwas_gwas AS g
  ON g.id = d.SNPID
LEFT JOIN mcpwas_pheno AS p
  ON p.study = d.study
  AND p.substudy = d.substudy
  AND p.phenotype = d.phenotype
LEFT JOIN mcpwas_class AS c
  ON c.phenotype = d.phenotype
LIMIT 10

