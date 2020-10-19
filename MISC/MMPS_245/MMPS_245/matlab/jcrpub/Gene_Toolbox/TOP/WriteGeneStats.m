function WriteGeneStats(filename,nom_pval,emp_pval,emp_adj_pval,snp_names,pheno_names)

numsnps = length(snp_names);
numpheno = length(pheno_names);

[fid,msg]=fopen(filename,'w');
if fid==-1
   error(msg);
end

fprintf(fid, 'SNP_name');


if 0

for ii=1:length(pheno_names)
   pstr = pheno_names{ii};
   fprintf(fid,', %s_nom_pval, %s_emp_pval, %s_emp_corr_pval, %s_nom_snp_x_sex, %s_emp_snp_x_sex, %s_emp_corr_snp_x_sex', pstr, pstr, pstr, pstr, pstr, pstr);
end
fprintf(fid,'\n');

for ii=1:numsnps,
   fprintf(fid,'%s',snp_names{ii});
   for jj=1:numpheno,
      fprintf(fid,', %.2e, %.2e, %.2e, %.2e, %.2e, %.2e', nom_pval(2,ii,jj), emp_pval(2,ii,jj), emp_adj_pval(2,ii,jj), nom_pval(8,ii,jj), emp_pval(8,ii,jj), emp_adj_pval(8,ii,jj));
   end
   fprintf(fid,'\n');
end

else

for ii=1:length(pheno_names)
   pstr = pheno_names{ii};
   fprintf(fid,', %s_all_nom, %s_all_emp, %s_all_emp_corr, %s_male_nom, %s_male_emp, %s_male_emp_corr, %s_female_nom, %s_female_emp, %s_female_emp_corr', pstr, pstr, pstr, pstr, pstr, pstr, pstr, pstr, pstr);
end
fprintf(fid,'\n');

for ii=1:numsnps,
   fprintf(fid,'%s',snp_names{ii});
   for jj=1:numpheno,
      fprintf(fid,', %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e', nom_pval(3,ii,jj), emp_pval(3,ii,jj), emp_adj_pval(3,ii,jj), nom_pval(2,ii,jj), emp_pval(2,ii,jj), emp_adj_pval(2,ii,jj), nom_pval(1,ii,jj), emp_pval(1,ii,jj), emp_adj_pval(1,ii,jj));
   end
   fprintf(fid,'\n');
end



end

fclose(fid);
