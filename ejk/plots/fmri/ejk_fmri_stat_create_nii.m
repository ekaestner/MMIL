% cd /home/mmilmcd/data/MCD_BOLD/subjects/EVENTfc028/orig.BLOCK/
% 3dmerge -dxyz=1 -1noneg -1clust 1 20 -1thresh 2.59 -1tindex 67 -1dindex 66 -prefix N-FF_cs20 stats.block.EVENTfc028+orig

function ejk_fmri_stat_create_nii(cfg)

if ~isfield(cfg,'pre_fix'); cfg.pre_fix = [ cfg.sbj_nme '_' cfg.fnc_mri_nme '_' 'Tstat' '_' cfg.mdl_fnc_nme ]; end

%% Find Volume
[~,nme_hld] = unix(['3dinfo -label' ' ' cfg.bld_dir '/' cfg.bld_fle_nme '.nii.gz']);
nme_hld = regexp(nme_hld,'\|','split'); nme_hld{end} = nme_hld{end}(1:end-1);

tst_nme_hld = find(strcmpi(nme_hld,[cfg.fnc_mri_stt '_Tstat']))-1;

%% Clean up prior stats


%% Make Stats
cmd = '';
cmd = [cmd 'cd' ' ' cfg.bld_dir ' ' ];
cmd = sprintf('%s\n',cmd);

bth_cmd     = [cmd '3dmerge' ' ' ...
                   '-dxyz=1' ' ' ...
                   '-1clust 1' ' ' num2str(cfg.cls_sze) ' ' ...
                   '-2thresh'  ' ' num2str(-cfg.sig_thr)  ' ' num2str(cfg.sig_thr) ' ' ...
                   '-1tindex'  ' ' num2str(tst_nme_hld) ' ' ...
                   '-1dindex'  ' ' num2str(tst_nme_hld) ' ' ...
                   '-prefix'   ' ' cfg.pre_fix '.nii.gz' ' ' ...
                   '-session'  ' ' cfg.out_dir ' ' ...
                   cfg.bld_fle_nme '.nii.gz' ];
out = unix(bth_cmd);

end

% pos_cmd     = [cmd '3dmerge' ' ' ...
%                    '-dxyz=1' ' ' ...
%                    '-1noneg' ' ' ...
%                    '-1clust 1' ' ' num2str(cfg.cls_sze) ' ' ...
%                    '-1thresh'  ' ' num2str(cfg.sig_thr) ' ' ...
%                    '-1tindex'  ' ' num2str(tst_nme_hld) ' ' ...
%                    '-1dindex'  ' ' num2str(tst_nme_hld) ' ' ...
%                    '-prefix'   ' ' cfg.sbj_nme '_' cfg.fnc_mri_nme '_' 'Tstat' '_' 'positive' '.nii.gz' ' ' ...
%                    '-session'  ' ' cfg.out_dir ' ' ...
%                    cfg.bld_fle_nme '.nii.gz' ];
% out = unix(pos_cmd);
%
%
% 