% 
% 
% 
% 

function ejk_fmri_surf_extract(fcfg)

%% Surf Out
ejk_chk_dir([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/'])

%
sbj_dir_lst = dir(sprintf('%s/FSURF*',fcfg.fsr_dir));
n = regexp({sbj_dir_lst.name},['FSURF_' fcfg.fsr_sbj_nme '.+_1$'],'match'); n = [n{:}];
fsr_hld = n{1};

%
if exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc = 'orig';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig.BLOCK' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc = 'orig.BLOCK';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig2' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc = 'orig2';
end

%
cmd = '';
cmd = sprintf('%scd %s\n',cmd,[fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/']);

if exist([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'lhs.1D'])==2; delete([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'lhs.1D']); end
if exist([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'rhs.1D'])==2; delete([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'rhs.1D']); end

lft_cmd     = [cmd '3dVol2Surf' ' ' ...
               '-spec ' fcfg.fsr_dir  '/' fsr_hld '/' 'SUMA' '/' fcfg.bld_sbj_nme '_both.spec' ' ' ...
               '-surf_A ' 'lh.' fcfg.f__mri_srf ' ' ...
               '-sv ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' fcfg.bld_sbj_nme '_SurfVol_Alnd_Exp+orig' ' ' ...
               '-grid_parent ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig' ' ' ... % '-grid_parent ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' '''stats.' fcfg.bld_sbj_nme '+orig[' num2str(nme_hld) ']''' ' ' ...
               '-map_func ' 'mask' ' ' ...
               '-out_1D ' fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'lhs.1D' ];
out = unix(lft_cmd);

rgh_cmd     = [cmd '3dVol2Surf' ' ' ...
               '-spec ' fcfg.fsr_dir  '/' fsr_hld '/' 'SUMA' '/' fcfg.bld_sbj_nme '_both.spec' ' ' ...
               '-surf_A ' 'rh.' fcfg.f__mri_srf ' ' ...
               '-sv ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' fcfg.bld_sbj_nme '_SurfVol_Alnd_Exp+orig' ' ' ...
               '-grid_parent ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig' ' ' ... % '-grid_parent ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' '''stats.' fcfg.bld_sbj_nme '+orig[' num2str(nme_hld) ']''' ' ' ...
               '-map_func ' 'mask' ' ' ...
               '-out_1D ' fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' fcfg.f__mri_srf '_' 'rhs.1D' ];
out = unix(rgh_cmd);

end