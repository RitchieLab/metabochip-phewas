load data local infile '~/group/projects/MetaboPheWAS/dbSNP/chr_1.txt'
into table _dbsnp_snp_import
ignore 7 lines
(rs,@_,@_,@_,@_,@_,chr,@_,@_,@_,@_,@pos,@_,@_,@_,@_,valid,@_,@_,@_,@_,source,@_,@_,@_,@_)
set pos=(case when @pos='' or @pos=' ' then null else @pos end), genes=''

load data local infile '~/group/projects/MetaboPheWAS/dbSNP/RsMergeArch.bcp'
into table _dbsnp_merge_import
(rs_old,rs_new,build,@_,@_,@_,rs_cur,@_,@_)

