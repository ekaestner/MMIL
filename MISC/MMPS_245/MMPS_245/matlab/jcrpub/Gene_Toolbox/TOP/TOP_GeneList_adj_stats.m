function [TOP_stats] = TOP_GeneList_adj_stats(GeneList,basedir,numperms,write_flag)


for ii = 1:length(GeneList),
      gene = GeneList{ii};
      gene_dir = sprintf('%s/%s',basedir,gene);
      nom_stats_file = sprintf('%s/stats/%s_nom_stats.mat',gene_dir,gene);
      outfile = sprintf('%s/stats/%s_adj_stats.mat',gene_dir,gene);
      if write_flag
         adj_stats_fname = outfile;
      else
         adj_stats_fname = '';
      end

      perm_str = sprintf('%dk',round(numperms/1000));
      null_stats_files.dir = sprintf('%s/stats',gene_dir);
      null_stats_files.name = sprintf('%s_null_stats_%s_snp*',gene,perm_str);

      [emp_pval,emp_adj_pval,pval_info] = TOP_gene_adj_stats(gene,nom_stats_file,null_stats_files,adj_stats_fname);

      TOP_stats(ii).gene_name = gene;
      TOP_stats(ii).info = pval_info;
      TOP_stats(ii).emp_pval = emp_pval;
      TOP_stats(ii).emp_adj_pval = emp_adj_pval;
end

