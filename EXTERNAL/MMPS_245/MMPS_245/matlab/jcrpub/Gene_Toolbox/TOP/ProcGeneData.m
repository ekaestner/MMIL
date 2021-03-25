function [stats] = ProcGeneData(gene,pheno_class,indatdir,projdir,sexchkfile,write_flag)

% ProcGeneData.m
%
% Usage: stats = ProcGeneData(gene,pheno_class,indatdir,projdir,write_flag)
%
% ProcGeneData returns regression statistics for association between
%  the specified phenotypes and GENE SNP's. The results are optionally
%  saved to a file if WRITE_FLAG is set.
%
% Inputs:
%   gene        - gene name (e.g. 'MECP2')
%   pheno_class - 'ROI' or 'COG'
%   indatdir    - input data dir name
%   projdir     - base project directory name
%   sexchkfile  - sexcheck results filename
%   write_flag  - force writing of nominal stats output file 
%
% Outputs:
%   stats       - struct containing regression statistics
%

snpfile = sprintf('%s/SNP_table_%s.csv',indatdir,gene);
if ~exist('snpfile','file')
  error(sprintf('Input file %s does not exist.',snpfile));
end

[SNP,Pheno,SubjInfo] = LoadGeneData(gene,pheno_class,indatdir,sexchkfile);

% Does the Proc<gene>_<pheno_class>_extra function exist for this gene? If so, run it.
gene_dir = sprintf('%s/%s/%s',projdir,pheno_class,gene);
custom_proc_fun = sprintf('Proc%s_%s_extra',gene,pheno_class);
if exist(sprintf('%s/%s.m',gene_dir,custom_proc_fun),'file'),
   cd(gene_dir);
   eval(sprintf('proc_further = @%s;',custom_proc_fun));
   [SNP,Pheno] = proc_further(SNP,Pheno);
end

% Load autosome_flag and regr struct for this gene
load(sprintf('%s/%s/%s/process_params.mat',projdir,pheno_class,gene));

% do some quality control on SNP and phenotype data
[SNP,Pheno,SubjInfo] = SNP_pheno_QC(SNP,Pheno,SubjInfo,autosome_flag);

BP_col = (SubjInfo.diagn==1);
SCZ_col = (SubjInfo.diagn==2);
%other_col = (SubjInfo.diagn > 2);
other_col = (SubjInfo.diagn==3 | SubjInfo.diagn==4);

males = (SubjInfo.sex==1);

CovData = [SubjInfo.sex SubjInfo.age BP_col SCZ_col other_col];
numvars = size(CovData,2)+1;

% Make x2xf matrix for regstats if needed (see x2xf.m for help)
if ischar(regr.model) && strcmp(regr.model,'x2fx_matrix')
   regr.model = [zeros(1,numvars); eye(numvars); [1 1 zeros(1,numvars-2)]];
end
%regr.group_sepvar = false; % snp term is for all subjects
%regr.group_sepvar = true;  % each sex will get its own snp term in design matrix
%regr.contrast = [0 1 1 0 0 0 0 0];  % mean of male and female x snp effects
%regr.tests = {'intercept', 'female_x_snp', 'male_x_snp', 'sex', 'age', 'BP', 'SCZ', 'other_diag', 'mean_male_female_x_snp'};
%regr.tests = {'intercept', 'snp', 'sex', 'age', 'BP', 'SCZ', 'other_diag', 'snp_x_sex'};

regr.groupID = males;
regr.numperms = 0;

% Perform regression specified in regr struct for all pair-wise combinations of columns of
% SNP.dat and Pheno.dat
stats = regstats_colwise(SNP.dat,CovData,Pheno.dat,regr);

% Only return stats for regressors that will undergo permutation test
if isempty(regr.contrast)
   keep_beta_ind = regr.perm_regressor_ind;
else
   keep_beta_ind = regr.perm_regressor_ind(1:end-1);
end
stats.beta = stats.beta(keep_beta_ind,:,:);
stats.beta_se = stats.beta_se(keep_beta_ind,:,:);
stats.t = stats.t(regr.perm_regressor_ind,:,:);
stats.pval = stats.pval(regr.perm_regressor_ind,:,:);
stats.tests = regr.tests(regr.perm_regressor_ind);

nom_stats_file = sprintf('%s/stats/%s_nom_stats.mat',gene_dir,gene);

if write_flag || ~exist(nom_stats_file,'file'), 
   cmd = sprintf('save %s -v7.3 SNP Pheno SubjInfo CovData regr stats',nom_stats_file);
   eval(cmd);
end
