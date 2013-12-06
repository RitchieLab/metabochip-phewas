DROP TABLE gwas_catalog;
CREATE TABLE gwas_catalog (
  rowid INTEGER PRIMARY KEY NOT NULL AUTO_INCREMENT,
  date_added DATE NOT NULL,
  pubmedID INT UNSIGNED NOT NULL,
  author VARCHAR(128) NOT NULL,
  date_published DATE NOT NULL,
  journal VARCHAR(64) NOT NULL,
  url VARCHAR(256) NOT NULL,
  study VARCHAR(512) NOT NULL,
  disease_trait VARCHAR(256) NOT NULL,
  initial_sample_size VARCHAR(1024) NOT NULL,
  replication_sample_size VARCHAR(1024) NOT NULL,
  region VARCHAR(64) NOT NULL,
  chr TINYINT UNSIGNED NULL,
  pos BIGINT UNSIGNED NULL,
  reported_genes VARCHAR(512) NOT NULL,
  mapped_genes VARCHAR(512) NOT NULL,
  upstream_entrezID INT UNSIGNED NULL,
  downstream_entrezID INT UNSIGNED NULL,
  SNP_entrezIDs VARCHAR(128) NOT NULL,
  upstream_gene_distance DOUBLE NULL,
  downstream_gene_distance DOUBLE NULL,
  SNP_risk_allele VARCHAR(64) NULL,
  SNPs VARCHAR(256) NULL,
  merged TINYINT NULL,
  SNP_current BIGINT UNSIGNED NULL,
  context VARCHAR(128) NOT NULL,
  intergenic TINYINT NULL,
  risk_AF VARCHAR(32) NOT NULL,
  pval DOUBLE NULL,
  pval_mlog DOUBLE NULL,
  pval_text VARCHAR(64) NOT NULL,
  ORbeta DOUBLE NULL,
  95_CI VARCHAR(64) NULL,
  platform VARCHAR(64) NULL,
  CNV VARCHAR(2) NOT NULL,
  KEY (chr,pos,rowid)
);

LOAD DATA LOCAL INFILE '/home/atf3/work/phewas/gwascatalog.txt'
INTO TABLE gwas_catalog
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(
  @date_added,pubmedID,author,@date_published,journal,url,
  study,disease_trait,initial_sample_size,@replication_sample_size,region,@chr,@pos,
  @reported_genes,mapped_genes,@upstream_entrezID,@downstream_entrezID,SNP_entrezIDs,@upstream_gene_distance,@downstream_gene_distance,
  @SNP_risk_allele,@SNPs,@merged,@SNP_current,context,@intergenic,@risk_AF,
  @pval,@pval_mlog,@pval_text,@ORbeta,@95_CI,platform,CNV
)
SET
  date_added=STR_TO_DATE(@date_added,'%m/%d/%Y'),
  date_published=STR_TO_DATE(@date_published,'%m/%d/%Y'),
  replication_sample_size=(CASE WHEN @replication_sample_size = 'NR' THEN '' ELSE @replication_sample_size END),
  chr=(CASE WHEN @chr = '' THEN NULL ELSE @chr END),
  pos=(CASE WHEN @pos = '' THEN NULL ELSE @pos END),
  reported_genes=(CASE WHEN @reported_genes = 'NR' THEN '' ELSE @reported_genes END),
  upstream_entrezID=(CASE WHEN @upstream_entrezID = '' THEN NULL ELSE @upstream_entrezID END),
  downstream_entrezID=(CASE WHEN @downstream_entrezID = '' THEN NULL ELSE @downstream_entrezID END),
  upstream_gene_distance=(CASE WHEN @upstream_gene_distance = '' THEN NULL ELSE @upstream_gene_distance END),
  downstream_gene_distance=(CASE WHEN @downstream_gene_distance = '' THEN NULL ELSE @downstream_gene_distance END),
  SNP_risk_allele=(CASE WHEN @SNP_risk_allele = 'NR' THEN '' ELSE @SNP_risk_allele END),
  SNPs=(CASE WHEN @SNPs = 'NR' THEN '' ELSE @SNPs END),
  merged=(CASE WHEN @merged = '' THEN NULL ELSE @merged END),
  SNP_current=(CASE WHEN @SNP_current = '' THEN NULL ELSE @SNP_current END),
  intergenic=(CASE WHEN @intergenic = '' THEN NULL ELSE @intergenic END),
  risk_AF=(CASE WHEN @risk_AF IN ('','NR') THEN '' ELSE @risk_AF END),
  pval=(CASE WHEN @pval IN ('','NS','E','Pending') THEN NULL WHEN @pval = '0E0' THEN 0.0 ELSE @pval END),
  pval_mlog=(CASE WHEN @pval_mlog IN ('','NS') THEN NULL ELSE @pval_mlog END),
  pval_text=(CASE WHEN @pval_text = 'NS' THEN '' ELSE @pval_text END),
  ORbeta=(CASE WHEN @ORbeta IN ('NR','#######','Pending') THEN NULL WHEN @ORbeta = '0E0' THEN 0.0 ELSE @ORbeta END),
  95_CI=(CASE WHEN @95_CI IN ('NR','[NR]') THEN '' WHEN @95_CI = '0E0' THEN 0.0 ELSE @95_CI END) 
;

SHOW WARNINGS;
