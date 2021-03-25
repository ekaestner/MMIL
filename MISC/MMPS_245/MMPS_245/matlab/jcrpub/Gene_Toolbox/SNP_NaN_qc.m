function [subj_qc_ok, SNP_qc_ok] = SNP_NaN_qc(SNP_mat)

% SNP_NaN_qc.m
% Usage: [subj_qc_ok, SNP_qc_ok] = SNP_NaN_qc(SNP_mat)
%
%  SNP_NaN_qc assumes that both the missing and invalid genotypes are
%  represented as NaN's in the SNP_mat matrix. SNP_mat is numsubj x numsnp.
%

numsubj = size(SNP_mat,1);
numsnp = size(SNP_mat,2);
snp_ind = [1:numsnp];
subj_ind = [1:numsubj];

% Find keeper SNPs/subjects
XX = isnan(SNP_mat);
keep_subj = [sum(XX,2) < ceil(0.5*numsnp)];
SNP_mat = SNP_mat(keep_subj,:);
subj_ind = subj_ind(keep_subj);
XX = isnan(SNP_mat);
keep_snp = find(sum(XX,1) < ceil(0.05*length(subj_ind)));
SNP_mat = SNP_mat(:,keep_snp);
%keep_subj = [sum(isnan(SNP_mat),2) == 0];
%SNP_mat = SNP_mat(keep_subj,:);

subj_qc_ok = false(numsubj,1);
subj_qc_ok(subj_ind) = true;

SNP_qc_ok = false(numsnp,1);
SNP_qc_ok(snp_ind(keep_snp)) = true;
