function [Te_ISA,Te,Tm_ISA,Tm] = trans_matrix_3shell(mode,test_opt,ISA,DIM,s,H_inv,P,A,no_megsensor)
%TRANS_MATRIX_3SHELL - calculate BEM transfer matrices for EEG and MEG
% function [Te_ISA,Te,Tm_ISA,Tm] = trans_matrix_3shell(mode,test_opt,ISA,DIM,s,H_inv,P,A,no_megsensor)
%**************************************************************************
%
% this function calculates transfer matrices for EEG and MEG
% mode        =1 , compute Te
%              2 , compute Tm
%              3 , compute BOTH
% test_opt    : 0, collocation,
%             : 1, Galerkin
% ISA         : 0 for NOT using isolated skull approach
%             : 1 for using the isolated skull approach
% DIM         : DIM(i) is the number of triangles(or nodes) on surface i
%             : 3 x 1
% s           : [conductivity of brain, conductivity of skull]
%             : 2 x 1
% P           : EEG scan matrix       no_eegsensor x sum(DIM)
% H_inv      : inverse of LU decomposition of EEG system matrix H
% A           : MEG system matrix
% no_megsensor: number of MEG sensors
%
% The isolated skull approach is incorporated for THREE shells model.
%*****************************************************************************
% last Modified by M.X. Huang, July 2006
% other people invloved in old version: J. Chang, J. Mosher, E. Ermer.
