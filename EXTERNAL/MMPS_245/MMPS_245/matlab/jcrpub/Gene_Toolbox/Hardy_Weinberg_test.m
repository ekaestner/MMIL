function [p, q, x2] = Hardy_Weinberg_test(geno_calls)

%  [p, q, x2] = Hardy_Weinberg_test(geno_calls)
%
% Only valid for autosomal and pseudo-autosomal SNP's

  n0 = sum(geno_calls==0);
  n1 = sum(geno_calls==1);
  n2 = sum(geno_calls==2);
  n = n0 + n1 + n2;
  p = (2*n0 + n1)/max(eps,2*n);
  q = (2*n2 + n1)/max(eps,2*n);
  n0_exp = n*p^2;
  n1_exp = 2*n*p*q;
  n2_exp = n*q^2;
  x2dev = (n0-n0_exp)^2/max(eps,n0_exp) + (n1-n1_exp)^2/max(eps,n1_exp) + (n2-n2_exp)^2/max(eps,n2_exp);
  x2 = 1-chi2cdf(x2dev,1);
