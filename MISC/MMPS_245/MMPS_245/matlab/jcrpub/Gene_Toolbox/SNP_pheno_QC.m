function [SNP,Pheno,SubjInfo] = SNP_pheno_QC(SNP,Pheno,SubjInfo,autosome_flag)

keep_snps = []; jj=0;
for ii=1:size(SNP.dat,2),

  if autosome_flag,
%    [p, qv, x2] = Hardy_Weinberg_test(SNP.dat(:,ii));
%    if p<.10 || p>0.90, continue; end
    maf = snp_MAF(SNP.dat(:,ii));
  else
    maf = chrX_snp_MAF(SNP.dat(:,ii), SubjInfo.sex);
  end
%  fprintf(1,'%s MAF = %.2f\n', SNP.snp_names{ii}, maf);
  %if maf*size(SNP.dat,1) < 30, continue; end
  if maf < 0.05, 
%     fprintf(1,'Deleting %s, MAF < 0.05\n', rsSNP_names{ii});
      continue;
  end

  jj = jj+1;
  keep_snps(jj) = ii;
end
SNP.dat = SNP.dat(:,keep_snps);
SNP.snp_names = SNP.snp_names(keep_snps);

[subj_qc_ok, SNP_qc_ok] = SNP_NaN_qc(SNP.dat);
SNP.dat = SNP.dat(subj_qc_ok, SNP_qc_ok);
SNP.snp_names = SNP.snp_names(SNP_qc_ok);
SNP.subjID = SNP.subjID(subj_qc_ok);
Pheno.subjID = Pheno.subjID(subj_qc_ok);
Pheno.dat = Pheno.dat(subj_qc_ok,:);
SubjInfo.ID = SubjInfo.ID(subj_qc_ok);
SubjInfo.sex = SubjInfo.sex(subj_qc_ok);
SubjInfo.age = SubjInfo.age(subj_qc_ok);
SubjInfo.eth = SubjInfo.eth(subj_qc_ok);
SubjInfo.diagn = SubjInfo.diagn(subj_qc_ok);

%if sum(isnan(SNP.dat(:))) > 0
%   error('Invalid SNPs still present after QC.');
%end

% Find keeper pheno's/subjects
[subj_qc_ok, pheno_qc_ok] = pheno_NaN_qc(Pheno.dat);
Pheno.dat = Pheno.dat(subj_qc_ok,pheno_qc_ok);
Pheno.subjID = Pheno.subjID(subj_qc_ok);
Pheno.names = Pheno.names(pheno_qc_ok);
SNP.dat = SNP.dat(subj_qc_ok,:);
SNP.subjID = SNP.subjID(subj_qc_ok);
SubjInfo.ID = SubjInfo.ID(subj_qc_ok);
SubjInfo.sex = SubjInfo.sex(subj_qc_ok);
SubjInfo.age = SubjInfo.age(subj_qc_ok);
SubjInfo.eth = SubjInfo.eth(subj_qc_ok);
SubjInfo.diagn = SubjInfo.diagn(subj_qc_ok);
