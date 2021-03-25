% cd /home/mmilmcd/data/MCD_BOLD/subjects/EVENTfc028/orig.BLOCK/
% 3dmerge -dxyz=1 -1noneg -1clust 1 20 -1thresh 2.59 -1tindex 67 -1dindex 66 -prefix N-FF_cs20 stats.block.EVENTfc028+orig

function ejk_fmri_stat_create(fcfg)

%% Setup
%
if     exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig' '/' 'stats.' fcfg.bld_sbj_nme '+orig.BRIK'])==2
    dta_loc = 'orig';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig.BLOCK' '/' 'stats.' fcfg.bld_sbj_nme '+orig.BRIK'])==2
    dta_loc = 'orig.BLOCK';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' 'orig2' '/' 'stats.' fcfg.bld_sbj_nme '+orig.BRIK'])==2
    dta_loc = 'orig2';
end

%
[~,nme_hld] = unix(['3dinfo -label' ' ' fcfg.bld_dir '/' fcfg.bld_sbj_nme '/' dta_loc '/' 'stats.' fcfg.bld_sbj_nme '+orig']);
nme_hld = regexp(nme_hld,'\|','split'); nme_hld{end} = nme_hld{end}(1:end-1);
nme_hld = find(strcmpi(nme_hld,fcfg.f__mri_stt))-1;

%% Create Stats
if exist(  [ fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    delete([ fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK']);
    delete([ fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.HEAD']);
end

cmd = '';
cmd = sprintf('%scd %s/%s/%s\n',cmd,fcfg.bld_dir,fcfg.bld_sbj_nme,dta_loc);

lft_cmd     = [cmd '3dmerge' ' ' ...
               '-dxyz=1' ' ' ...
               '-1noneg' ' ' ...
               '-1clust 1' ' ' num2str(fcfg.cls_sze) ' ' ...
               '-1thresh' ' ' num2str(fcfg.sig_thr) ' ' ...
               '-1tindex' ' ' num2str(nme_hld) ' ' ...
               '-1dindex' ' ' num2str(nme_hld) ' ' ...
               '-prefix' ' ' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig' ' ' ...
               '-session ' [ fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' ] ' ' ...
               'stats.' fcfg.bld_sbj_nme '+orig.BRIK' ];
out = unix(lft_cmd);

cmd = '';
cmd = sprintf('%scd %s\n',cmd,[ fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'surf_Out' ]);
[~,out] = unix([cmd 'gunzip' ' ' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK.gz']);

end