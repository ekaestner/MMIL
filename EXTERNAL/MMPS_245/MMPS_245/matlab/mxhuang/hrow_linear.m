function Hi = hrow_linear(colp,PP,v,yset,yset_mag,geometry,nodes,cdvdiff,N,area2)
%HROW_LINEAR - Compute row in BEM geometry matrix for linear collocation
% function Hi = hrow_linear(colp,PP,v,yset,yset_mag,geometry,nodes,cdvdiff,N,area2)
% compute a row of the H matrix(1 by no_node) for linear BEM
% formula in:
% de Munck, IEEE Trans.BME, pp986-990,1992
% or schlitt et al, IEEE Trans, BME, pp52-58, 1995
% some global variables
%
% last modofocation: M.X. Huang, MArch 2004
% Other people invloved in old version: J. Chang, J. Mosher, E. Ermer

