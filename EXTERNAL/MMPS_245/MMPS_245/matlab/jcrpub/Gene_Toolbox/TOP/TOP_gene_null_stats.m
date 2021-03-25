%OlesGeneList = {'EPHA4','EPHB2','DISC1','CACNG2','RELN','PDE4B','MDGA1','NTRK3','GRIK3','GRIK4','MCTP2','DCLK1','BRCA2','TNFRSF1A','RIPK3','NMB','LPA','B3GAT2','HAPLN2','IL1F8','FAS','OLEAND1'};

%GeneList = {'MECP2','MCPH1','CDK5RAP2','ASPM','CENPJ',}; 

%GeneList = {'PDE4B','MDGA1','NTRK3','GRIK3','GRIK4','MCTP2','DCLK1','BRCA2','TNFRSF1A','RIPK3','NMB','LPA','B3GAT2','HAPLN2','IL1F8','FAS','OLEAND1'};

%GeneList = {'CDK5RAP2','ASPM','CENPJ','EPHA4','EPHB2','DISC1','CACNG2','RELN','PDE4B','MDGA1','NTRK3','GRIK3','GRIK4','MCTP2','DCLK1','BRCA2','TNFRSF1A','RIPK3','NMB','LPA','B3GAT2','HAPLN2','IL1F8','FAS','OLEAND1','MCPH1'};

cluster_startup;

GeneList = {'ASPM','CENPJ','CDK5RAP2','MCPH1'};

write_flag = 1;
snps_per_outfile = 3;
numperms = 100000;
pheno_class = 'COG';  % 'COG' or 'ROI'

TOPdir = '/home/cooper/space_md8_5/TOPgeno_md8_5';
TOP_stats_dir = '/home/cooper/space_md8_5/TOPgeno_md8_5/GeneStats';
TOP_doc_dir = '~dale/TOP_docs';

sexchkfile = '/home/cooper/space_md8_5/TOPgeno_md8_5/TOP_GW_6_0.sexcheck';

for ii = 1:length(GeneList),

   gene = GeneList{ii};
   gene_dir = sprintf('%s/%s/%s',TOP_stats_dir,pheno_class,gene);

   if write_flag || ~exist(nom_stats_file,'file')
      ProcGeneData(gene,pheno_class,TOP_doc_dir,TOP_stats_dir,sexchkfile,write_flag);   % creates nom_stats_file
   end

   nom_stats_file = sprintf('%s/stats/%s_nom_stats.mat',gene_dir,gene);
   load(nom_stats_file);  % load stats struct
   num_snps = size(stats.pval,2);

   cmd = sprintf('%s/Gene_q_null_segs.csh %s %s %d %d %d %s',TOPdir,TOPdir,gene,num_snps,snps_per_outfile,numperms,gene_dir);

   system(cmd);
end


% The following can only be run after all submitted cluster jobs finish.
% Convert nominal association stats to empirical stats and stats corrected for multiple comparisons.
if 0
   basedir = sprintf('%s/%s',TOP_stats_dir,pheno_class);
   [TOP_stats] = TOP_GeneList_adj_stats(GeneList,basedir,numperms,write_flag);
end 
