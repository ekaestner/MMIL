function [Rstats] = regstats_colwise(X,CovData,Y,regr)

% regstats_colwise.m
%
% restats_colwise() performs linear regression on all possible pairwise combi-
% nations of the columns of input matrix X and output matrix Y. CovData is the
% covariate matrix included in each regression. The numbers of rows in each of
% these matrices must agree.
%
% Inputs:
%   X       - independent variable matrix (numsamples x numvariables)
%   CovData - covariate data matrix (numsamples x numcovars)
%   Y       - dependent variable matrix ((numsamples x numdepvars)
%   regr    - struct specifying details of the regression
%
% Output:
%   Rstats  - struct containing t-stats, p-vals, and beta coeffs (only for
%             non-permuted case.
% 
% regr struct:
% regr.model    - design matrix specification. See help for regstats.m.
% regr.contrast - contrast vector (length must agree with number of columns in
%                 the regressor matrices that will be created from X and CovData).
% regr.numperms - if exists and > 0, number of permutations to perform
% regr.groupID  - vector of group ID's. Length must equal the number of rows of X.
% regr.group_sepvar - if true, prior to running regressions, separate each column
%                     of X into multiple regressor columns, w/ one column per
%                     unique group value in regr.groupID.
%
% restats_colwise can be used to generate distributions of null-hypothesis stats
%   via permutations of the rows of X. Empirical estimates of significance can be
%   made from nominal stats vs. the distribution of null-hypothesis stats.
%
% For a regression model involving a given response (Y_col_j) and predictor matrix
% [X_col_i CovData], if a NaN is present in a row of the predictor matrix or the 
% corresponding element of Y_col_j, then that row and element will be removed.
%

if ~isequal(size(X,1),size(Y,1))
   error('X and Y have unequal number of rows.');
end

if isfield(regr,'numperms') && regr.numperms > 0
   perm_flag = true;
   numiters = regr.numperms;
else
   regr.numperms = 0;
   perm_flag = false;
   numiters = 1;
end

stats_list = {'tstat', 'mse', 'covb'};

if ~isfield(regr,'contrast')
   regr.contrast = [];
end

if ~isfield(regr,'model')
   error('REGR struct arg must contain "model" field.');
end

if  ~isfield(regr,'group_sepvar')
   regr.group_sepvar = false;
end
if regr.group_sepvar & ~isfield(regr,'groupID')
   error('If regr.group_sepvar is true, then regr must have "groupID" field also.');
end

if isfield(regr,'groupID')
   groups = unique(regr.groupID);
   for ii=1:length(groups),
      group_ind(ii) = {find(regr.groupID==groups(ii))};
   end
else
   % effectively create a single 'group' by including all samples in one
   group_ind = {[1:size(X,1)]};
end
numgroups = length(group_ind);

num_indep_vars = 1;
if regr.group_sepvar, num_indep_vars = numgroups; end

% which stats from regstats should be saved?
if perm_flag,
   % Only regstats() output for the indpendent vars will be saved.
   % The 1st var for regstats is the intercept, which we skip.
   keepterms = 1+[1:num_indep_vars];
else
   % No permutation, so keep all nominal stats and alloc mem for beta coeffs.
   numterms = 1+num_indep_vars+size(CovData,2);
   keepterms = [1:numterms];
   Rstats.beta = zeros(numterms,size(X,2),size(Y,2));
   Rstats.beta_se = zeros(size(Rstats.beta));
end

nkterms = length(keepterms);  % terms to save from regstats()
test_terms = keepterms;       % keepterms and possibly a constrast term
if ~isempty(regr.contrast), 
   test_terms = [test_terms length(regr.tests)]; 
end

Rstats.tests = regr.tests(test_terms);

% stats for contrast (if present) will be included in Rstats.t and Rstats.pval
Rstats.t = zeros(length(test_terms),size(X,2),size(Y,2),numiters);
Rstats.pval = zeros(size(Rstats.t));

orig_ind_seq = [1:size(X,1)];

for np=1:numiters,

   new_inds = orig_ind_seq;  % no permutation

   if perm_flag
      if ~rem(np,100), fprintf(1,'%d  ',np); end

      if regr.group_sepvar,
         new_inds = zeros(size(X,1),1);

         % permute samples within each group, not across groups
         for ng=1:numgroups,
            tt = group_ind{ng};
            new_inds(tt) = tt(randperm(length(tt)));
         end

      else
         % standard permutation
         new_inds = randperm(size(X,1));
      end
   end

   % Loop over columns of X
   for ii=1:size(X,2),
      xvec = X(new_inds,ii);  % copy of possibly permuted column of X
      Xmat = zeros(size(X,1),numgroups);  % independent variable(s) 
      for ng=1:numgroups,
         inds = group_ind{ng};
         Xmat(inds,ng) = xvec(inds,1);
      end

      Cmat = [Xmat CovData];   % tack on unpermuted covariate columns

      % Loop over columns of Y (dependent vars), performing the specified regression
      for jj=1:size(Y,2),
         D = regstats(Y(:,jj),Cmat,regr.model,stats_list);
         Rstats.t(1:nkterms,ii,jj,np) = D.tstat.t(keepterms);
         Rstats.pval(1:nkterms,ii,jj,np) = D.tstat.pval(keepterms);
         if ~isempty(regr.contrast)
            xtxi = D.covb/D.mse;
            sigma_contrast = sqrt(regr.contrast*xtxi*regr.contrast');
            contrast_dotp = regr.contrast*D.tstat.beta;
            contrast_t = contrast_dotp/(sigma_contrast*sqrt(D.mse));
            Rstats.t(end,ii,jj,np) = contrast_t;
            Rstats.pval(end,ii,jj,np) = 2*(tcdf(-abs(contrast_t), D.tstat.dfe));
         end
    
         % Save beta coeff and std errs if this is nominal stats
         if ~perm_flag,
           Rstats.beta(:,ii,jj) = D.tstat.beta;
           Rstats.beta_se(:,ii,jj) = D.tstat.se;
         end
      end
   end
end
fprintf(1,'\n');

