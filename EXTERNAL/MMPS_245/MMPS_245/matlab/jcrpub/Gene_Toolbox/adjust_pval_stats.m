function [emp_pval,emp_adj_pval] = adjust_pval_stats(pval,null_logpval)

% null_logpval MUST be num_predictors x num_permutations

if ~isvector(pval) || length(pval) ~= size(null_logpval,1)
   fprintf(1,'Sizes of input args are either invalid or mismatched');
   return;
end
numperms = size(null_logpval,2);
if numperms < 1000,
   error('This analysis probably will fail for small number of permutations.')
end

% For each nominal p-val (pval), first determine the corresponding 'empirical'
% p-val by scaling the nom. p-val by the slope of the best fit line (passing 
% through the origin) of the cdf-derived permuted p-vals vs. the permuted
% p-vals (null_logpval). The null-hypoth distrib for the ith nominal p-val
% is approximated by hist(null_logpval(i,:))/size(null_logpval,2).
% This procedure should coerce the scaled nominal p-vals to have approx. the
% same standard dev. Next, correct each empirical p-val for the multiple
% comparisons ('experiment-wise' correction) by first determining the 
% distrib of the min (across comparisons) of the slope-scaled null-hypoth
% p-vals. For each empirical p-val, then determine the proportion of scaled
% null-hypoth p-vals that are lower. This proportion, emp_adj_pval(i), is
% the estimate of the true, corrected p-value for the ith predictor of the
% dependent variable.

adj_slope = zeros(size(pval));
emp_pval = zeros(size(pval));
scaled_logp = zeros(size(null_logpval));

for ii=1:length(pval),

%   null_logp = squeeze(null_logpval(ii,:));
%   standard_emp_pval(ii) = length(find(null_logp > -log10(pval(ii))))/numperms;

   logp = squeeze(null_logpval(ii,:));
   [tvhist,bctr] = hist(logp,numperms);
   tvcdf = cumsum(tvhist)/numperms;
   bctr_range = max(bctr)-min(bctr);
   lowbval = min(bctr) + 0.2*bctr_range;
   hibval = max(bctr) - 0.25*bctr_range;
   min_ind = find(bctr>lowbval,1);
   max_ind = find(bctr>=hibval,1);
   if isempty(min_ind) || isempty(max_ind)
      error('Empty xmin_ind or xmax_ind.');
   end
   x = bctr(min_ind:max_ind)';
   y = -log10(1-tvcdf(min_ind:max_ind))';
   adj_slope(ii) = x\y;
   scaled_logp(ii,:) = adj_slope(ii)*logp';
   emp_pval(ii) = 10^-(-log10(pval(ii))*adj_slope(ii));
end

max_logp = max(scaled_logp,[],1);  % max across the numsnps dim
emp_adj_pval = zeros(size(pval));
for ii=1:length(emp_pval),
   emp_adj_pval(ii) = length(find(max_logp > -log10(emp_pval(ii))))/numperms;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
for ii=1:numsnps,
   tstr = pstats.snp_names{ii};
   tstr(strfind(tstr,'"')) = '';
   pstats.snp_names(ii) = {tstr};
end
for jj=1:numpheno,
   tstr = pstats.pheno_names{jj};
   tstr(strfind(tstr,'"')) = '';
   pstats.pheno_names(jj) = {tstr};
end

fid=fopen('/home/cooper/TOPgeno/MECP2/Alex_snps/MECP2_pval_table_10k_1.csv','w');
fprintf(fid, 'SNP_name'); 
for ii=1:length(pstats.pheno_names)
   pstr = pstats.pheno_names{ii};
   fprintf(fid,', %s_nom_pval_male, %s_empir_pval_male, %s_empir_adj_pval_male, %s_nom_pval_female, %s_empir_pval_female, %s_empir_adj_pval_female', pstr, pstr, pstr, pstr, pstr, pstr);
end
fprintf(fid,'\n');
for ii=1:numsnps,
   fprintf(fid,'%s',pstats.snp_names{ii});
   for jj=1:numpheno,
      fprintf(fid,', %.2e, %.2e, %.2e, %.2e, %.2e, %.2e', stats.pval(beta_ind(1),ii,jj), emp_pval(1,ii,jj), emp_adj_pval(1,ii,jj), stats.pval(beta_ind(2),ii,jj), emp_pval(2,ii,jj), emp_adj_pval(2,ii,jj));
   end
   fprintf(fid,'\n');
end
fclose(fid);


fid=fopen('/home/cooper/TOPgeno/MECP2/MECP2_pval_table.csv','w');
fprintf(fid, 'SNP_name'); fprintf(fid,', %s', pstats.pheno_names{:});
for ii=1:length(pstats.pheno_names)
   pstr = pstats.pheno_names{ii};
   fprintf(fid,', %s_nom_pval_all, %s_empir_pval_all, %s_empir_adj_pval_all, %s_nom_pval_male, %s_empir_pval_male, %s_empir_adj_pval_male, %s_nom_pval_female, %s_empir_pval_female, %s_empir_adj_pval_female', pstr, pstr, pstr, pstr, pstr, pstr, pstr, pstr, pstr);
end
fprintf(fid,'\n');
for ii=1:numsnps,
   fprintf(fid,'%s',pstats.snp_names{ii});
   for jj=1:numpheno,
      fprintf(fid,', %.2e, %.2e, %.2e', stats.pval_all(ii,jj), stats.emp_pval_all(ii,jj), stats.emp_adj_pval_all(ii,jj), stats.pval_male(ii,jj), stats.emp_pval_male(ii,jj), stats.emp_adj_pval_male(ii,jj), stats.pval_female(ii,jj), stats.emp_pval_female(ii,jj), stats.emp_adj_pval_female(ii,jj));
   end
   fprintf(fid,'\n');
end
fclose(fid);
end
