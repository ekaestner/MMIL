function ejk_fishersZ( cfg )

%% Save Data
for iC = 1:size(cfg.rvl_one,2)
    rvl_one.(['rvl_one' num2str(iC)]) = cfg.rvl_one(:,iC);
    rvl_two.(['rvl_two' num2str(iC)]) = cfg.rvl_two(:,iC);
    num_one.(['num_one' num2str(iC)]) = cfg.num_one(:,iC);
    num_two.(['num_two' num2str(iC)]) = cfg.num_two(:,iC);
end

ejk_chk_dir(cfg.out_dir)

save( [ cfg.out_dir '/' 'rvl_one.mat' ], 'rvl_one' )
save( [ cfg.out_dir '/' 'rvl_two.mat' ], 'rvl_two' )
save( [ cfg.out_dir '/' 'num_one.mat' ], 'num_one' )
save( [ cfg.out_dir '/' 'num_two.mat' ], 'num_two' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_fishersZ.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_fishersZ( ''%s/rvl_one.mat'', ''%s/rvl_two.mat'',''%s/num_one.mat'', ''%s/num_two.mat'', ''%s'', ''%s'')', sve_cmd, ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), cfg.out_pre );
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' ejk_fix_path(cfg.out_dir) '/example_R_script.r' ] );

%% Tidy up
fsh_out = mmil_readtext( [ cfg.out_dir '/' cfg.out_pre '_fishersZ.csv' ] );
fsh_out{1,1} = '';
fsh_out(1,2:end) = cfg.col_lbl;
fsh_out(2:end,1) = cfg.row_lbl;
cell2csv( [ cfg.out_dir '/' cfg.out_pre '_fishersZ.csv' ], fsh_out );

zou_out = mmil_readtext( [ cfg.out_dir '/' cfg.out_pre '_Zou2007.csv' ] );
zou_out{1,1} = '';
zou_out(1,2:end) = cfg.col_lbl;
zou_out(2:end,1) = cfg.row_lbl;
cell2csv( [ cfg.out_dir '/' cfg.out_pre '_Zou2007.csv' ], zou_out );


end