function pvl_out = ejk_surface_pvalue_cluster( cfg , pvl_inp )

tmp_dir = '/home/ekaestne/Downloads';

tmp_fle_inp = [ tmp_dir '/' 'pvl_thr_input.mgh' ];
tmp_fle_out = [ tmp_dir '/' 'pvl_thr_output.w' ];

%%
pvl_inv = 1-pvl_inp;
thr_inv = 1-cfg.pvl_thr;

fs_save_mgh( pvl_inv , tmp_fle_inp );

%%
cmd = 'mri_surfcluster';
cmd = sprintf('%s --in %s'               ,cmd, tmp_fle_inp);
cmd = sprintf('%s --thmin %0.6f'         ,cmd, thr_inv);
cmd = sprintf('%s --subject %s'          ,cmd, cfg.fsr_sbj);
cmd = sprintf('%s --hemi %s'             ,cmd, cfg.hms);
cmd = sprintf('%s --sd %s'               ,cmd, cfg.fsr_dir);
cmd = sprintf('%s --o %s'                ,cmd, tmp_fle_out);
cmd = sprintf('%s --sum %s'              ,cmd, [ tmp_fle_out(1:end-2) '.txt' ]);
cmd = sprintf('%s --minarea %0.6f'       ,cmd, cfg.cls_thr);
[status,result] = unix(cmd);

%%
[ msk_val , msk_ind ] = fs_read_wfile(tmp_fle_out);
msk_ind = msk_ind(find(msk_val)); % exclude 0 vals (shouldn't be there anyway)

pvl_out          = zeros(size(pvl_inv));
pvl_out(msk_ind) = pvl_inv(msk_ind);

pvl_out = 1 - pvl_out;

%%
delete(tmp_fle_inp);
delete(tmp_fle_out);
  
end