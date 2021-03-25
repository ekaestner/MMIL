TOPdir = '~/space_md8_5/TOPgeno_md8_5';
%outdir = sprintf('%s/TOP_docs',TOPdir);
outdir = sprintf('%s/TOP_docs/Lars_mcph_csv',TOPdir);

%GeneList = {'EPHA4','EPHB2','DISC1','CACNG2','RELN','PDE4B','MDGA1','NTRK3','GRIK3','GRIK4','MCTP2','DCLK1','BRCA2','TNFRSF1A','RIPK3','NMB','LPA','B3GAT2','HAPLN2','IL1F8','FAS','OLEAND1'};

%GeneList = {'ASPM','CENPJ','CDK5RAP2','MCPH1'};
GeneList = {'ASPM','CENPJ','CDK5RAP2','MCPH1'};
pheno_class = '';  % 'COG' or 'ROI'

for ii=1:length(GeneList)

   gene = GeneList{ii};
   fprintf(1,'%s ',gene);
   matfile=sprintf('%s/GeneStats/%s/%s/stats/%s_adj_stats.mat',TOPdir,pheno_class,gene,gene);
   load(matfile);
   outfile = sprintf('%s/%s_%s_stats.csv',outdir,pheno_class,gene);
   WriteGeneStats(outfile,stats.pval,emp_pval,emp_adj_pval,SNP.snp_names,Pheno.names);
end
fprintf(1,'\n');
