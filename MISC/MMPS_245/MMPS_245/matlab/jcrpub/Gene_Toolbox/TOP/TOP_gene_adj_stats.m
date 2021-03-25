function [emp_pval,emp_adj_pval,pval_info] = TOP_gene_adj_stats(gene,nom_stats_file,null_stats_files,outfile)

%  TOP_gene_adj_stats.m
%
%  Usage: [emp_pval,emp_adj_pval,pval_info] = TOP_gene_adj_stats(gene,nom_stats_file,null_stats_files,outfile)
%
%  TOP_gene_adj_stats converts nominal p-values for SNP's in GENE to empirical p-values
%   and p-values corrected for multiple comparisons.
%
%  Inputs:
%    gene             - gene name
%    nom_stats_file   - file containing nominal p-values
%    null_stats_files - struct has the fields: .dir (directory name) and .name (a regular expression)
%    outfile          - output filename (if empty, stats not saved to file)
%
%  Outputs:
%    emp_pval         -  empirical p-values
%    emp_adj_pval     -  p-values corrected for multiple comparisons
%    pval_info        -  struct containing tests, snp, and phenotype names
%

% load nominal stats
if ~exist(nom_stats_file,'file')
   error(sprintf('File %s not found.',nom_stats_file))
end
load(nom_stats_file);

total_snps = size(SNP.dat,2);
numphenos = size(Pheno.dat,2);
numterms = size(stats.pval,1);

% Force stats.pval to be 3-D
stats.pval = reshape(stats.pval, [numterms total_snps numphenos]);

% Get names of all our null stats files
D=dir(sprintf('%s/%s',null_stats_files.dir,null_stats_files.name));
if isempty(D)
   error(sprintf('No null stats files found for %s.',gene));
end
cd(null_stats_files.dir);

% Load first null stats file and get necessary info
tmp=load(D(1).name);
sz_logpval = size(tmp.pstats.logpval);
nperms =  sz_logpval(end);
regr.numperms = nperms;
logpval=zeros([size(stats.pval) nperms],'single');

% Load all null log(pvals)
ind1=1;
for ii=1:length(D),
   tmp=load(D(ii).name);
   fprintf(1,'%s loaded.\n', D(ii).name);

   % Confirm pstats.logpval is 4-D
   if ~isequal(ndims(tmp.pstats.logpval),4),
      error('pstats.logpval is not 4-D.');
%     tmp.pstats.logpval = single(reshape(tmp.pstats.logpval,[numterms num_snp numphenos nperms]));
   end

   num_snp = size(tmp.pstats.logpval,2); % num. of snp's in the current null_stats file

   ind2=ind1-1+num_snp;

   logpval(:,ind1:ind2,:,:) = tmp.pstats.logpval;
  
   ind1=ind1+num_snp;
end

emp_pval = zeros(length(stats.tests),size(SNP.dat,2),size(Pheno.dat,2));
emp_adj_pval = zeros(length(stats.tests),size(SNP.dat,2),size(Pheno.dat,2));
fprintf(1,'%s: ',gene);
for nr=1:length(stats.tests),
   fprintf(1,'%d ',nr);
   for np=1:size(Pheno.dat,2),
      null_logp = squeeze(logpval(nr,:,np,:));

      % adjust_pval_stats() expects that the 1st dim of null_logp is the snp dim. If
      %total_snps=1, then the above squeeze() call introduces a transpose which must
      % be reversed.
      if total_snps==1,  null_logp=null_logp';  end
     
      [emp_pval(nr,:,np),emp_adj_pval(nr,:,np)] = adjust_pval_stats(colvec(stats.pval(nr,:,np)),null_logp);
   end
end
fprintf(1,'\n');

% So we know what each p-value in emp_pval and emp_adj_pval means...
pval_info.tests = stats.tests;          % 1st dimension
pval_info.snp_names = SNP.snp_names;    % 2nd dim.
pval_info.pheno_names = Pheno.names;    % 3rd dim.

if ~isempty(outfile),
   cmd = sprintf('save %s -v7.3 SNP Pheno SubjInfo CovData emp_adj_pval emp_pval stats regr',outfile);
   eval(cmd);
end
