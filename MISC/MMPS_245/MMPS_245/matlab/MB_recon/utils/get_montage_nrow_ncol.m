function [nrow, ncol] = get_montage_nrow_ncol(n)
%
% function [nrow, ncol] = get_montage_nrow_ncol(n)
%
% Set up number of columns and rows for displaying n images, using code adapted from montage.m
%
% Input
%   n    - Number of images
%
% Outputs
%   nrow - Number of rows in montage
%   ncol - Number of columns in montage
%
% (c) Kangrong Zhu      Stanford University     June 2014

aspectRatio = 1; % 1 for a roughly square display
ncol = sqrt(aspectRatio * n);
ncol = ceil(ncol);
nrow = ceil(n / ncol);

return