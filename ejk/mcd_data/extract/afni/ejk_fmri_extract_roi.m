% fcfg = [];
% 
% fcfg.prj_dir = prj_dir; '/home/ekaestne/PROJECTS/'
% fcfg.prj_nme = prj_nme; 'fMRISurfImages';
% 
% fcfg.sbj_nme = sbj_nme{iS};
% 
% fcfg.fsr_dir = fsr_dir; '/home/mmilmcd/data/FSRECONS/';
% fcfg.bld_dir = bld_dir; '/home/mmilmcd/data/MCD_BOLD/subjects/';
% 
% fcfg.fsr_sbj_nme = fsr_sbj_nme; 
% fcfg.bld_sbj_nme = bld_sbj_nme; 
% 
% fcfg.f__mri_stt = f__mri_stt{iT}; { 'N-FF_GLT#0_Tstat'  'N#0_Tstat'         'FF#0_Tstat' };
% fcfg.f__mri_nme = f__mri_nme{iT}; { 'N_FF'              'NW'                'FF' };
% 
% fcfg.f__mri_typ = f_mri_typ{iFT};
% 
% fcfg.lbl_nme = lbl_nme{iL}; {'' '.a2009s'}
% 
% fcfg.cls_sze = cls_sze; 20;

function ejk_fmri_extract_roi(fcfg)

%% Setup
if exist([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'Background' '/' 'ROI_Holder']) ~= 7; mkdir([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'Background' '/' 'ROI_Holder']); end

fsr_col = fcfg.fsr_sbj_nme{1,2}; fsr_col = regexp(fsr_col,';','split'); fsr_col = find(strcmpi(fsr_col,'fMRI'));
fsr_sbj_row = find(strcmpi(fcfg.fsr_sbj_nme(:,1),fcfg.sbj_nme));
fsr_hld = fcfg.fsr_sbj_nme{fsr_sbj_row,2}; fsr_hld = regexp(fsr_hld,';','split'); fsr_hld = fsr_hld{fsr_col};

roi_nme_hld = mmil_readtext([fcfg.fsr_dir '/' fsr_hld '/' 'label' '/' 'aparc.annot' fcfg.lbl_nme '.ctab'],'  ');
roi_nme_hld = [roi_nme_hld(1:10,2:3) ; roi_nme_hld(11:end,1:2)];
lhs_roi_nme_hld = strcat('ctx_lh_',roi_nme_hld(:,2));
rhs_roi_nme_hld = strcat('ctx_rh_',roi_nme_hld(:,2));

bld_sbj_row = find(strcmpi(fcfg.bld_sbj_nme(:,1),fcfg.sbj_nme));
if exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig' '/' 'aparc' fcfg.lbl_nme '+aseg_rank_Alnd_Exp+orig.BRIK'])==2
    dta_loc = 'orig';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig.BLOCK' '/' 'aparc' fcfg.lbl_nme '+aseg_rank_Alnd_Exp+orig.BRIK'])==2
    dta_loc = 'orig.BLOCK';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig2' '/' 'aparc' fcfg.lbl_nme '+aseg_rank_Alnd_Exp+orig.BRIK'])==2
    dta_loc = 'orig2';
else
     error(['Subject Error : ' fcfg.sbj_nme])
end

if exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc2 = 'orig';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig.BLOCK' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc2 = 'orig.BLOCK';
elseif exist([fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' 'orig2' '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.BRIK'])==2
    dta_loc2 = 'orig2';
else
     error(['Subject Error : ' fcfg.sbj_nme])
end

%% Pull ROI    
cmd = '';
cmd     = sprintf('%scd %s/%s/%s\n',cmd,fcfg.bld_dir,fcfg.bld_sbj_nme{bld_sbj_row,2},dta_loc);
cmd     = [cmd '3dROIstats' ' ' ...
           '-nomeanout' ' ' ...
           fcfg.f__mri_typ ' ' ...
           '-mask' ' ' ...
           'aparc' fcfg.lbl_nme '+aseg_rank_Alnd_Exp+orig.' ' ' ...
           fcfg.bld_dir '/' fcfg.bld_sbj_nme{bld_sbj_row,2} '/' dta_loc2 '/' fcfg.f__mri_nme '_' 'cs' num2str(num2str(fcfg.cls_sze)) '+orig.'];
[~,out] = unix(cmd);
    
out_reg = regexp(out,'\s','split'); out_reg = out_reg(1:end-1);
top_row = out_reg(1:numel(out_reg)/2); top_row = strrep(top_row,'NZcount_','');
bot_row = out_reg(numel(out_reg)/2+1:end);
out_reg = [top_row' bot_row'];

out_reg(:,1) = strrep(out_reg(:,1),'-','_');

[~,out_reg_lhs_int,lhs_int] = intersect(out_reg(:,1),lhs_roi_nme_hld);
[~,out_reg_rhs_int,rhs_int] = intersect(out_reg(:,1),rhs_roi_nme_hld);

lhs_out = [lhs_roi_nme_hld num2cell(zeros(numel(lhs_roi_nme_hld),1))];
lhs_out(lhs_int,2) = num2cell(cellfun(@str2num,out_reg(out_reg_lhs_int,2)));

rhs_out = [rhs_roi_nme_hld num2cell(zeros(numel(rhs_roi_nme_hld),1))];
rhs_out(rhs_int,2) = num2cell(cellfun(@str2num,out_reg(out_reg_rhs_int,2)));

out_fle = [lhs_out ;  rhs_out];

cell2csv([fcfg.prj_dir '/' 'DATA' '/' fcfg.sbj_nme '/' 'BOLD' '/' 'ROI'  '/' fcfg.sbj_nme '_' fcfg.f__mri_nme '_' 'aparc' mmil_spec_char(fcfg.lbl_nme,{'.'}) '_' mmil_spec_char(fcfg.f__mri_typ,{'-'}) '.csv'],out_fle)

end
